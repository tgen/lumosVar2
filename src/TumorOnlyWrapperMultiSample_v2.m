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
%delete(gcp('nocreate'));
%distcomp.feature( 'LocalUseMpiexec', true);
%parpool(inputParam.numCPU);
% pc=parcluster();
% pc.NumWorkers=inputParam.numCPU;
% saveAsProfile(pc,'pc');
% parallel.defaultClusterProfile('pc');

Tcell=cell('');
Ecell=cell('');
%%%% If outmat exits, uses outmat, otherwise generates data tables
if(exist([inputParam.outMat],'file'))
    vars={'T','E','Tcell','Ecell'};
    load([inputParam.outMat],'-mat',vars{:});
    inputParam=readInputs(paramFile)
else
   [T, E]=preprocessTumorOnly_v2(inputParam,paramFile);
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

pool=gcp('nocreate');
if isempty(pool)
    delete(gcp('nocreate'));
    %distcomp.feature( 'LocalUseMpiexec', true);
    pc = parcluster('local');
    pc.NumWorkers = inputParam.numCPU;
    parpool(pc, pc.NumWorkers);
end

for i=1:size(Ecell,2)
    E=Ecell{i};
    exonRD{i}=E{median(MapQC,2)<inputParam.minExonQual & median(perReadPass,2)>inputParam.minPerReadPASS & median(abFrac,2)>inputParam.minABFrac & min(TumorRD,[],2)>=0 & min(NormalRD,[],2)>=0,:};
    segs{i}=segmentData(exonRD{i},inputParam.cnaAlpha);
    clear E;
end
segsMerged=mergeSegments_v2(segs,exonRD,inputParam);

%%%Make sure segments extend to ends of chromosome
for i=1:22
    idx1=find(segsMerged(:,1)==i,1,'first');
    idx2=find(segsMerged(:,1)==i,1,'last');
    segsMerged(idx1,2)=min([Tcell{1}.Pos(Tcell{1}.Chr==i); exonRD{1}(exonRD{1}(:,1)==i,2)]);
    segsMerged(idx2,3)=max([Tcell{1}.Pos(Tcell{1}.Chr==i); exonRD{1}(exonRD{1}(:,1)==i,3)]);
end


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



%%%Fit Copy Number Model
pool=gcp('nocreate');
if isempty(pool)
    delete(gcp('nocreate'));
    %distcomp.feature( 'LocalUseMpiexec', true);
    pc = parcluster('local');
    pc.NumWorkers = inputParam.numCPU;
    parpool(pc, pc.NumWorkers);
end
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
    idx=getPosInRegionSplit([T.Chr T.Pos],segsTable{i}(:,1:3),inputParam.blockSize);
    T.NumCopies=segsTable{i}(idx,5);
    T.MinAlCopies=segsTable{i}(idx,6);
    T.cnaF=segsTable{i}(idx,7);
    T.W=segsTable{i}(idx,8);
    T.BmeanBQ(T.BCountF+T.BCountR==0)=inputParam.defaultBQ;
    Tcell{i}=T;
    clear T;
end
save([inputParam.outMat]);

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
    Tcell{j}.ApopAFcomb=countsAll.ApopAF(locb);
    Tcell{j}.BpopAFcomb=countsAll.BpopAF(locb);
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
save([inputParam.outMat]);
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
while(sum(abs(somPos-somPosOld))/sum(somPos)>0.05 && i<inputParam.maxIter)
     [segsTable, W, f, CNAscale, nll, t{i}]=fitCNAmulti(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam);
    %['clonal fractions: ' num2str(f)]
    for j=1:length(Tcell)
        T=Tcell{j};
        idx=getPosInRegionSplit([T.Chr T.Pos],segsTable{:,1:3},inputParam.blockSize);
        T.NumCopies=segsTable.N(idx);
        T.MinAlCopies=segsTable.M(idx);
        T.cnaF=segsTable.F(idx,j);
        Tcell{j}=T;
    end
    save([inputParam.outMat]);
    i=i+1
    somPosOld=somPos;
    inputParam.numClones=size(f,2);
    if inputParam.NormalSample>0
        n=1;
        NormalSample=inputParam.NormalSample;
        inputParam.NormalSample=1;
        inputParam.numClones=1;
        for j=1:length(Tcell)
            if j~=NormalSample
                idx=[NormalSample j];
                %postComb=jointSNV_v2(Tcell(idx), exonRD(idx), f(n,:), W(idx,:), inputParam);
                [postCombPair{n},~,pDataCombPair{n}]=jointSNV_v2(Tcell(idx), exonRD(idx), inputParam.priorF, 3*ones(2,1), inputParam);
                [lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[postCombPair{n}.Chr postCombPair{n}.Pos],'rows');
                P.SomaticPair(:,j)=postCombPair{n}.Somatic(locb);
                n=n+1;
            else
                P.SomaticPair(:,j)=zeros(size(P,1),1);
            end
        end
        inputParam.numClones=size(f,2);
        inputParam.NormalSample=NormalSample;
    else
        P.SomaticPair=zeros(size(P,1),length(Tcell));
    end
    [postComb, pDataSum(i+1),pDataComb,clones,prior,countsAll]=jointSNV_v2(Tcell, exonRD, f, W, inputParam);
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
        T=Tcell{j};
        E=Ecell{j};
        homPos=P.Hom(:,j)>0.5  & P.SomaticPair(:,j)<0.5 | (max([P.Somatic P.SomaticPair],[],2)>0.5 & (max(P.DataSomatic(:,j),P.DataNonDip(:,j))<=P.DataHom(:,j) | (T.BcountsComb==0 & T.A==T.RefComb))) | (P.Somatic(:,j)>0.5 & inputParam.NormalSample==j);
        P.homPos(:,j)=homPos;
        [F{j},postTrust,postArtifact]=qualDiscrimCalls(T,E,homPos,inputParam);
        P.trust(:,j)=postTrust(:,2);
        P.artifact(:,j)=postArtifact(:,1);
        clear T E;
    end
    filtPos=max(P.trust,[],2)>inputParam.pGoodThresh & max(P.artifact,[],2)<inputParam.pGoodThresh;
    hetPos=max(P.Het,[],2)>inputParam.pGermlineThresh & filtPos;
    somPos=max([P.Somatic P.SomaticPair],[],2)>inputParam.pSomaticThresh & filtPos & min([Tcell{1}.ApopAFcomb Tcell{1}.BpopAFcomb],[],2)<inputParam.maxSomPopFreq;
    if inputParam.NormalSample>0
        for j=1:tIdx
            somPos(max(P.SomaticPair(:,tIdx(j)),[],2)>inputParam.pSomaticThresh & min(P.trust(:,[tIdx(j) inputParam.NormalSample]),[],2)>inputParam.pGoodThresh & min([Tcell{1}.ApopAF Tcell{1}.BpopAF],[],2)<inputParam.maxSomPopFreq  &  max(P.artifact(:,[tIdx(j) inputParam.NormalSample]),[],2)<inputParam.pGoodThresh,:)=1;
        end
    end
    ['Somatic positions: ' num2str(sum(somPos))]
    ['Het positions: ' num2str(sum(hetPos))] 
end

['converged at iteration' num2str(i)]

%%%write output files

% 
% fAll=zeros(length(Tcell),inputParam.numClones);
% 
% fAll(tIdx,1:end)=reshape(f,[],inputParam.numClones);


cloneCounts=hist(cloneId(somPos,1),1:size(f,2));
[~,ord]=sort(cloneCounts,'descend');
fSort=f(:,ord);
[~,cloneIdSort]=ismember(cloneId,ord);
ord=[ord size(f,2)+1];
[~,segsTable.cnaIdx]=ismember(segsTable.cnaIdx,ord);

save([inputParam.outMat]);
[Filter,somaticDetected]=writeJointVCF(Tcell,P,fSort,cloneIdSort,F,inputParam);
writeSegVCF(segsTable,exonRD,CNAscale,Tcell,hetPos,inputParam)
writeCloneSummary(segsTable,Ecell,Tcell,fSort,cloneIdSort,inputParam,Filter,somaticDetected);
plotTumorOnly(exonRD,segsTable,CNAscale,fSort,Tcell,somPos,hetPos,cloneIdSort,inputParam)

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
