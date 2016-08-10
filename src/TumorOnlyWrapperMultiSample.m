function output=TumorOnlyWrapperMultiSample(paramFile,varargin)
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
    exonRD{i}=E{E.MapQC<inputParam.minExonQual & E.perReadPass>inputParam.minPerReadPASS & E.abFrac>inputParam.minABFrac & E.TumorRD>=0 & E.NormalRD>=0,:};
    segs{i}=segmentData(exonRD{i},inputParam.cnaAlpha);
    clear E;
end
segsMerged=mergeSegments(segs,Ecell,inputParam);


save([inputParam.outMat]);

%%%Quality Filtering
for i=1:size(Tcell,2)
    T=Tcell{i};
    E=Ecell{i};
    P{i}=table();
    [F{i},P{i}.postTrust,P{i}.postArtifact]=qualDiscrim(T,E,inputParam)
    filtPos{i}=P{i}.postTrust(:,2)>inputParam.pGoodThresh & P{i}.postArtifact(:,1)<inputParam.pGoodThresh;
    clear T E;
    ['PASS positions: ' num2str(sum(filtPos{i}))]
end
    
%%%Preliminary Variant Classification
for i=1:size(Tcell,2)
    T=Tcell{i};
    hetPos=min([T.ApopAF T.BpopAF],[],2)>inputParam.minHetPopFreq & filtPos{i} & T.BCountF+T.BCountR>=inputParam.minBCount;
    somPos=T.CosmicCount>1 & min([T.ApopAF T.BpopAF],[],2)<inputParam.maxSomPopFreq & filtPos{i};
    dataHet{i}=[T.Chr(hetPos) T.Pos(hetPos) T.ControlRD(hetPos) T.ReadDepthPass(hetPos) T.BCountF(hetPos)+T.BCountR(hetPos)];
    dataSom{i}=[T.Chr(somPos) T.Pos(somPos) T.ControlRD(somPos) T.ReadDepthPass(somPos) T.BCountF(somPos)+T.BCountR(somPos)];
    ['Somatic positions: ' num2str(sum(somPos))];
    ['Het positions: ' num2str(sum(hetPos))];
    clear T;
end

%%%Make sure segments extend to ends of chromosome
for j=1:size(Tcell,2)
    T=Tcell{j};
    for i=1:22
        idx1=find(segs{j}(:,1)==i,1,'first');
        idx2=find(segs{j}(:,1)==i,1,'last');
        segs{j}(idx1,2)=min([T.Pos(T.Chr==i); exonRD{j}(exonRD{j}(:,1)==i,2)]);
        segs{j}(idx2,3)=max([T.Pos(T.Chr==i); exonRD{j}(exonRD{j}(:,1)==i,3)]);
    end
    clear T;
end

%%%Fit Copy Number Model
for i=1:size(Tcell,2)
    [segsTable{i}, W{i}, f{i}, c{i}, nll{i}]=fitCNA(dataHet{i},dataSom{i},exonRD{i},segs{i},inputParam);
    ['clonal fractions: ' num2str(f{i})]
end

%%%Add copy number info to data table
for i=1:size(Tcell,2)
    T=Tcell{i};
    %idxExon=getPosInRegionSplit([T.Chr T.Pos],exonRD(:,1:3),inputParam.blockSize);
    idx=getPosInRegions([T.Chr T.Pos],segsTable{i}(:,1:3));
    T.NumCopies=segsTable{i}(idx,5,1);
    T.MinAlCopies=segsTable{i}(idx,6,1);
    T.cnaF=segsTable{i}(idx,7,1);
    T.W=segsTable{i}(idx,8,1);
    T.BmeanBQ(T.BCountF+T.BCountR==0)=inputParam.defaultBQ;
    Tcell{i}=T;
    clear T;
end

%%%Initial Bayesian Variant Calling
pDataSum(1)=0;
[postComb, pDataSum(2),somaticFlag,pDataComb,clones]=jointSNV(Tcell, exonRD, f(1,:), W(1,:), inputParam);
for j=1:length(Tcell)
    [lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[postComb.Chr postComb.Pos],'rows');
    P{1,j}.Somatic=postComb.Somatic(locb).*somaticFlag(locb,j);
    P{1,j}.Het=postComb.Het(locb);
    P{1,j}.Hom=postComb.Hom(locb);
    P{1,j}.DataSomatic=pDataComb{j}.Somatic(locb);
    P{1,j}.DataHet=pDataComb{j}.Het(locb);
    P{1,j}.DataHom=pDataComb{j}.Hom(locb);
    cloneId{1,j}=clones(locb);
end

%%%repeat fitting and variant calling until converges
i=1;
while(round(pDataSum(i+1))~=round(pDataSum(i)) && i<inputParam.maxIter)
    for j=1:length(Tcell)
        T=Tcell{j};
        hetPos=P{i,j}.Het>inputParam.pGermlineThresh & filtPos{j};
        somPos=P{i,j}.Somatic>inputParam.pSomaticThresh & filtPos{j} & min([T.ApopAF T.BpopAF],[],2)<inputParam.maxSomPopFreq;
        dataHet=[T.Chr(hetPos) T.Pos(hetPos) T.ControlRD(hetPos) T.ReadDepthPass(hetPos) T.BCountF(hetPos)+T.BCountR(hetPos)];
        dataSom=[T.Chr(somPos) T.Pos(somPos) T.ControlRD(somPos) T.ReadDepthPass(somPos) T.BCountF(somPos)+T.BCountR(somPos)];
        ['Somatic positions: ' num2str(sum(somPos))]
        ['Het positions: ' num2str(sum(hetPos))]
        [segsTable{i+1,j}, W{i+1,j}, f{i+1,j}, c{i+1,j}, nll{i+1,j}]=fitCNA(dataHet,dataSom,exonRD{j},segs{j},inputParam);
        ['clonal fractions: ' num2str(f{i+1,j})]
        idx=getPosInRegions([T.Chr T.Pos],segsTable{i+1,j}(:,1:3));
        T.NumCopies=segsTable{i+1,j}(idx,5);
        T.MinAlCopies=segsTable{i+1,j}(idx,6);
        T.cnaF=segsTable{i+1,j}(idx,7);
        T.W=segsTable{i+1,j}(idx,8);
    end
    %save([inputParam.outMat]);
    i=i+1
    [postComb, pDataSum(i+1),somaticFlag,pDataComb,clones]=jointSNV(Tcell, exonRD, f(i,:), W(i,:), inputParam);  
    for j=1:length(Tcell)
        [lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[postComb.Chr postComb.Pos],'rows');
        P{i,j}.Somatic=postComb.Somatic(locb).*somaticFlag(locb,j);
        P{i,j}.Het=postComb.Het(locb);
        P{i,j}.Hom=postComb.Hom(locb);
        P{i,j}.DataSomatic=pDataComb{j}.Somatic(locb);
        P{i,j}.DataHet=pDataComb{j}.Het(locb);
        P{i,j}.DataHom=pDataComb{j}.Hom(locb);
        cloneId{i,j}=clones(locb,j);
    end  
end

['converged at iteration' num2str(i)]

%%%write output files


save([inputParam.outMat]);
for i=1:length(Tcell)
    writeVCF(Tcell{i},P{end,i},f{end,i},W{end,i},P{1,i}.postTrust(:,2),P{1,i}.postArtifact(:,2),cloneId{end,i},F{i},inputParam,i);
end
    
%writeSegVCF(segsTable(:,:,end),inputParam);
%message='finished writing VCFs'
%writeCloneSummary(segsTable(:,:,end),E,T,pSomatic(:,end),postTrust(:,2),postArtifact(:,2),f(end,:),W(end,:),cloneId(:,end),inputParam);
%message='finished writing summary table'
%plotTumorOnly(exonRD,segsTable(:,:,end),c(end),f(end,:),T,somPos,hetPos,cloneId(:,end),inputParam);

close all;
output=0;
exit;
