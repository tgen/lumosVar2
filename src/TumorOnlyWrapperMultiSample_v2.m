function output=TumorOnlyWrapperMultiSample_v2(paramFile,varargin)
%TumorOnlyWrapper - Entry function for tumor only caller
%uses a bayesian framework to call somatic variants and germline variants
%from tumor only exome sequencing
%
% Syntax:  output = TumorOnlyWrapper(paramFile)
%
% Inputs:
%    paramFile - parameter file in yaml format, see configTemplate.yaml 
%
% Outputs:
%    output - returns 0 upon completion
%
% Example: 
%   TumorOnlyWrapper('sampleConfig.yaml')
%
% Other m-files required: callCNA.m, callSNV.m, fitCNA.m, nllCNA.m, 
%   plotTumorOnly.m, preprocessTumorOnly.m, segmentData.m, 
%   writeCloneSummary.m  writeSegVCF.m  writeVCF.m
% Other requirements: parsePileupData.packed.pl, samtools, htslib
% Subfunctions: none
% MAT-files required: cghcbshybridnu.mat

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

inputParam=readInputs(paramFile);
cd(inputParam.workingDirectory);
addpath(genpath(inputParam.workingDirectory));

%%% start parallel pool
delete(gcp('nocreate'));
distcomp.feature( 'LocalUseMpiexec', true);
parpool(inputParam.numCPU);

Tcell=cell('');
Ecell=cell('');
%%%% If outmat exits, uses outmat, otherwise generates data tables
if(exist([inputParam.outMat],'file'))
    vars={'T','E','Tcell','Ecell'};
    load([inputParam.outMat],'-mat',vars{:});
    inputParam=readInputs(paramFile)
else
   [T, E]=preprocessTumorOnly(inputParam,paramFile);
   save([inputParam.outMat]);
end

if isempty(Tcell)
    Tcell=T;
end
if isempty(Ecell)
    Ecell=E;
end


%%% Filters Exon Data and Segments
for i=1:size(Ecell,2)
    E=Ecell{i};
    MapQC(:,i)=E.MapQC;
    perReadPass(:,i)=E.perReadPass;
    abFrac(:,i)=E.abFrac;
    TumorRD(:,i)=E.TumorRD;
    NormalRD(:,i)=E.NormalRD;
end

for i=1:size(Ecell,2)
    E=Ecell{i};
    exonRD{i}=E{median(MapQC,2)<inputParam.minExonQual & median(perReadPass,2)>inputParam.minPerReadPASS & median(abFrac,2)>inputParam.minABFrac & min(TumorRD,[],2)>=0 & min(NormalRD,[],2)>=0,:};
    segs{i}=segmentData(exonRD{i},inputParam.cnaAlpha);
    clear E;
end
segsMerged=mergeSegments(segs,Ecell,inputParam);


save([inputParam.outMat]);

%%%Quality Filtering
P=table();
for i=1:size(Tcell,2)
    T=Tcell{i};
    E=Ecell{i};
    idx1=T.ApopAF+T.BpopAF>1 & T.ApopAF<T.BpopAF;
    T.BpopAF(idx1)=1-T.ApopAF(idx1)-3*inputParam.pvFreqIndel;
    idx2=T.ApopAF+T.BpopAF>1 & T.ApopAF>=T.BpopAF;
    T.ApopAF(idx2)=1-T.BpopAF(idx2)-3*inputParam.pvFreqIndel;
    Tcell{i}=T;
    [F{i},postTrust,postArtifact]=qualDiscrim(T,E,inputParam);
    P.trust(:,i)=postTrust(:,2);
    P.artifact(:,i)=postArtifact(:,1);
    clear T E;
    %['PASS positions: ' num2str(sum(filtPos{i}))]
end
filtPos=max(P.trust,[],2)>inputParam.pGoodThresh & max(P.artifact,[],2)<inputParam.pGoodThresh;
    
%%%Preliminary Variant Classification
if inputParam.NormalSample>0
    hetPos=min([Tcell{inputParam.NormalSample}.ApopAF Tcell{inputParam.NormalSample}.BpopAF],[],2)>inputParam.minHetPopFreq & filtPos & Tcell{inputParam.NormalSample}.BCountF+Tcell{inputParam.NormalSample}.BCountR>=inputParam.minBCount;
else
    for i=1:length(Tcell)
        Bcounts(:,i)=Tcell{i}.BCountF+Tcell{i}.BCountR;
    end
    hetPos=min([Tcell{1}.ApopAF Tcell{1}.BpopAF],[],2)>inputParam.minHetPopFreq & filtPos & max(Bcounts,[],2)>=inputParam.minBCount;
end    
somPos=Tcell{1}.CosmicCount>1 & min([Tcell{1}.ApopAF Tcell{1}.BpopAF],[],2)<inputParam.maxSomPopFreq & filtPos;

%%%Make sure segments extend to ends of chromosome
for i=1:22
    idx1=find(segsMerged(:,1)==i,1,'first');
    idx2=find(segsMerged(:,1)==i,1,'last');
    segsMerged(idx1,2)=min([Tcell{1}.Pos(Tcell{1}.Chr==i); exonRD{1}(exonRD{1}(:,1)==i,2)]);
    segsMerged(idx2,3)=max([Tcell{1}.Pos(Tcell{1}.Chr==i); exonRD{1}(exonRD{1}(:,1)==i,3)]);
end


%%%Fit Copy Number Model
inputParam.numClones=1;
parfor i=1:size(Tcell,2)
    dataHet=[Tcell{i}.Chr(hetPos) Tcell{i}.Pos(hetPos) Tcell{i}.ControlRD(hetPos) Tcell{i}.ReadDepthPass(hetPos) Tcell{i}.BCountF(hetPos)+Tcell{i}.BCountR(hetPos)];
    if i==inputParam.NormalSample
        dataSom=[];
    else
        dataSom=[Tcell{i}.Chr(somPos) Tcell{i}.Pos(somPos) Tcell{i}.ControlRD(somPos) Tcell{i}.ReadDepthPass(somPos) Tcell{i}.BCountF(somPos)+Tcell{i}.BCountR(somPos)];
    end
    [segsTable{i}, W{i}, f{i}, c{i}, nll{i}, ~]=fitCNA(dataHet,dataSom,exonRD{i},segsMerged,inputParam);
    ['clonal fractions ' num2str(i) ': ' num2str(f{i})]
end
%[segsTable{1}, W{1}, f{1}, c{1}, nll{1}]=fitCNAmulti(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam);


%%%Add copy number info to data table
for i=1:size(Tcell,2)
    T=Tcell{i};
    %idxExon=getPosInRegionSplit([T.Chr T.Pos],exonRD(:,1:3),inputParam.blockSize);
    idx=getPosInRegions([T.Chr T.Pos],segsTable{i}(:,1:3));
    T.NumCopies=segsTable{i}(idx,5);
    T.MinAlCopies=segsTable{i}(idx,6);
    T.cnaF=segsTable{i}(idx,7);
    T.W=segsTable{i}(idx,8);
    T.BmeanBQ(T.BCountF+T.BCountR==0)=inputParam.defaultBQ;
    Tcell{i}=T;
    clear T;
end

%%%Initial Bayesian Variant Calling
pDataSum(1)=0;
[postComb, pDataSum(2),somaticFlag,pDataComb,clones,prior,countsAll]=jointSNVinit(Tcell, exonRD, cell2mat(f)', cell2mat(W)', inputParam);
for j=1:length(Tcell)
    [lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[postComb.Chr postComb.Pos],'rows');
    Tcell{j}.RefComb=countsAll.Ref(locb);
    Tcell{j}.Acomb=countsAll.A(locb);
    Tcell{j}.Bcomb=countsAll.B(locb);
    Tcell{j}.AcountsComb=countsAll.Acounts(locb,j);
    Tcell{j}.BcountsComb=countsAll.Bcounts(locb,j);
    P.Somatic(:,j)=postComb.Somatic(locb).*somaticFlag(locb,j);
    P.Het(:,j)=postComb.Het(locb);
    P.Hom(:,j)=postComb.Hom(locb);
    P.NonDip(:,j)=postComb.NonDip(locb);
    P.DataSomatic(:,j)=pDataComb{j}.Somatic(locb);
    P.DataHet(:,j)=pDataComb{j}.Het(locb);
    P.DataHom(:,j)=pDataComb{j}.Hom(locb);
    P.DataNonDip(:,j)=pDataComb{j}.NonDip(locb);
    P.priorSomatic(:,j)=prior.Somatic(locb);
    P.priorHet(:,j)=prior.Het(locb);
    P.priorHom(:,j)=prior.Hom(locb);
    P.priorNonDip(:,j)=prior.nonDiploid(locb);
    cloneId(:,j)=clones(locb);
end


for i=1:length(Tcell)
    T=Tcell{i};
    E=Ecell{i};
    homPos=P.Hom(:,i)>0.5 | (P.Somatic(:,i)>0.5 & (max(P.DataSomatic(:,i),P.DataNonDip(:,i))<=P.DataHom(:,i) | (T.BcountsComb==0 & T.A==T.RefComb)));
    [F{i},postTrust,postArtifact]=qualDiscrimCalls(T,E,homPos,inputParam);
    P.trust(:,i)=postTrust(:,2);
    P.artifact(:,i)=postArtifact(:,1);
    clear T E;
end
filtPos=max(P.trust,[],2)>inputParam.pGoodThresh & max(P.artifact,[],2)<inputParam.pGoodThresh;


%%%repeat fitting and variant calling until converges
i=1;
somPosOld=somPos;
hetPos=max(P.Het,[],2)>inputParam.pGermlineThresh & filtPos;
somPos=max(P.Somatic,[],2)>inputParam.pSomaticThresh & filtPos & min([Tcell{1}.ApopAF Tcell{1}.BpopAF],[],2)<inputParam.maxSomPopFreq;
['Somatic positions: ' num2str(sum(somPos))]
['Het positions: ' num2str(sum(hetPos))]
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
while(sum(abs(somPos-somPosOld))>0 && i<inputParam.maxIter)
     [segsTable, W, f, CNAscale, nll, t{i}]=fitCNAmulti(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam);
    %['clonal fractions: ' num2str(f)]
    for j=1:length(Tcell)
        T=Tcell{j};
        idx=getPosInRegions([T.Chr T.Pos],segsTable{:,1:3});
        T.NumCopies=segsTable.N(idx);
        T.MinAlCopies=segsTable.M(idx);
        T.cnaF=segsTable.F(idx,j);
        Tcell{j}=T;
    end
    save([inputParam.outMat]);
    i=i+1
    inputParam.numClones=size(f,2);
    [postComb, pDataSum(i+1),pDataComb,clones,prior,countsAll]=jointSNV(Tcell, exonRD, f, W, inputParam);
    if inputParam.NormalSample>0
        n=1;
        NormalSample=inputParam.NormalSample;
        inputParam.NormalSample=1;
        for j=1:length(Tcell)
            if j~=NormalSample
                idx=[NormalSample j];
                postComb=jointSNV(Tcell(idx), exonRD(idx), f(n,:), W, inputParam);
                [lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[postComb.Chr postComb.Pos],'rows');
                P.SomaticPair(:,j)=postComb.Somatic(locb);
                n=n+1;
            end
        end
        inputParam.NormalSample=NormalSample;
    else
        P.SomaticPair=zeros(size(P,1),length(Tcell));
    end
    for j=1:length(Tcell)
        [lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[postComb.Chr postComb.Pos],'rows');
        Tcell{j}.RefComb=countsAll.Ref(locb);
        Tcell{j}.Acomb=countsAll.A(locb);
        Tcell{j}.Bcomb=countsAll.B(locb);
        Tcell{j}.AcountsComb=countsAll.Acounts(locb,j);
        Tcell{j}.BcountsComb=countsAll.Bcounts(locb,j);
        P.Somatic(:,j)=postComb.Somatic(locb);
        P.Het(:,j)=postComb.Het(locb);
        P.Hom(:,j)=postComb.Hom(locb);
        P.NonDip(:,j)=postComb.NonDip(locb);
        P.DataSomatic(:,j)=pDataComb{j}.Somatic(locb);
        P.DataHet(:,j)=pDataComb{j}.Het(locb);
        P.DataHom(:,j)=pDataComb{j}.Hom(locb);
        P.DataNonDip(:,j)=pDataComb{j}.NonDip(locb);
        P.priorSomatic(:,j)=prior.Somatic(locb);
        P.priorHet(:,j)=prior.Het(locb);
        P.priorHom(:,j)=prior.Hom(locb);
        P.priorNonDip(:,j)=prior.nonDiploid(locb);
        cloneId(:,j)=clones(locb);
    end
    somPosOld=somPos;
    hetPos=max(P.Het,[],2)>inputParam.pGermlineThresh & filtPos;
    somPos=max([P.Somatic P.SomaticPair],[],2)>inputParam.pSomaticThresh & filtPos & min([Tcell{1}.ApopAF Tcell{1}.BpopAF],[],2)<inputParam.maxSomPopFreq;
    ['Somatic positions: ' num2str(sum(somPos))]
    ['Het positions: ' num2str(sum(hetPos))] 
end

['converged at iteration' num2str(i)]

%%%write output files
fAll=zeros(length(Tcell),inputParam.numClones);

fAll(tIdx,1:end)=reshape(f,[],inputParam.numClones);

save([inputParam.outMat]);
writeJointVCF(Tcell,P,f,cloneId,F,inputParam)
writeSegVCF(segsTable,inputParam)
writeCloneSummary(segsTable,Ecell,Tcell,P,f,cloneId,inputParam)
plotTumorOnly(exonRD,segsTable,CNAscale,f,Tcell,somPos,hetPos,cloneId,inputParam)

%for i=1:length(Tcell)
%writeVCF(Tcell{i},P,fAll(i,:),cloneId{i},F{i},inputParam,i);
% end
    
%writeSegVCF(segsTable(:,:,end),inputParam);
%message='finished writing VCFs'
%writeCloneSummary(segsTable(:,:,end),E,T,pSomatic(:,end),postTrust(:,2),postArtifact(:,2),f(end,:),W(end,:),cloneId(:,end),inputParam);
%message='finished writing summary table'
%plotTumorOnly(exonRD,segsTable(:,:,end),c(end),f(end,:),T,somPos,hetPos,cloneId(:,end),inputParam);

close all;
output=0;
exit;
