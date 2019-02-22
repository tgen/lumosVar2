function output=lumosVarMain(paramFile,numCPU)
%/Volumes/home/rhalperin/lumosVarTallArray/lumosVar2/src/lumosVarMain.m
%TumorOnlyWrapper - Entry function for tumor only caller
%uses a bayesian framework to call somatic variants and germline variants
%from tumor only exome sequencing
%
% Syntax:  output = TumorOnlyWrapper(paramFile)
%
% Inputs:
%    paramFile - parameter file in yaml format, see configTemplate.yaml
%    numCPU - number of CPU's to use for multithreading (pass as string)
%
% Outputs:
%    output - returns 0 upon completion
%
% Example: 
%   LumosVarMain('sampleConfig.yaml','16')
%
% Other m-files required: callCNAmulti.m, callVariants.m, jointSNV.m, 
%   fitCNAmulit.m, nllCNAmulti.m, nllCNAaddClone.m 
%   plotTumorOnly.m, preprocessTumorOnly.m, segmentData.m, 
%   writeCloneSummary.m  writeSegVCF.m  writeJointVCF.m, readInputs.m
% Other requirements: bcftools
% Subfunctions: none
% MAT-files required: cghcbshybridnu.mat

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-May-2018

%------------- BEGIN CODE --------------

['paramFile is ' paramFile{:}] 
profile -memory on;
setpref('profiler','showJitLines',1);

addpath(pwd())
inputParam=readInputs(paramFile);
cd(inputParam.workingDirectory);

%%%% Initial parallel pool

if nargin==2
    inputParam.numCPU=str2num(numCPU);
end
pool=gcp('nocreate');
%if isempty(pool)
    delete(gcp('nocreate'));
    pc = parcluster('local');
    pc.NumWorkers = inputParam.numCPU;
    parpool(pc, pc.NumWorkers, 'IdleTimeout', Inf);
%end

Tcell=cell('');
Ecell=cell('');
%%%% If outmat exits, uses outmat, otherwise generates data tables
if(exist([inputParam.outMat],'file'))
    vars={'posTable','E','Tcell','Ecell','sampleNames','bamList','inputParam'};
    load([inputParam.outMat],'-mat',vars{:});
    NormalSample=inputParam.NormalSample;
    priorF=inputParam.priorF;
    chrList=inputParam.chrList;
    sampleCount=inputParam.sampleCount;
    if exist('sampleNames')==0
        sampleNames=strsplit(inputParam.sampleNames,',');
    end
    inputParam=readInputs(paramFile);
    inputParam.NormalSample=NormalSample;
    inputParam.chrList=chrList;
    inputParam.sampleCount=sampleCount;
    if length(priorF)==length(Tcell)
        inputParam.priorF=priorF;
    end
else
    [posTable, Ecell,inputParam]=readBams(inputParam,paramFile{:});   
    fid=fopen(inputParam.bamList);
    bamList=textscan(fid,'%s');
    sampleCount=length(bamList{1});
    fclose(fid);
    for i=1:sampleCount
        paths=regexp(bamList{1}{i},',','split');
        folders=regexp(paths{1},'/','split');
        names=regexprep(folders{end},'\.bam$','');
        sampleNames{i}=names;
        message=['sample ' num2str(i) ': ' sampleNames{i} ' from ' bamList{1}{i}]
    end
    p=profile('info');
    %writetable(p.FunctionTable,[inputParam.outName 'prof.txt']);
    profile resume -history
    save([inputParam.outMat],'-v7.3');
end


%%%% get chromosomes to examine

%%%% find positions within bed file that have a population allele frequency
%%%% greater than maxSomPopFreq
parfor i=1:length(inputParam.chrList)
    cmd=strcat({'bcftools view -q '},num2str(inputParam.maxSomPopFreq),{' -R '},inputParam.regionsFile,{' '},inputParam.snpVCFpath,inputParam.chrList(i),inputParam.snpVCFname,' | grep -v ^# | awk ''{ print $2 "\n" }''')
    [~,output]=system(cmd{1});
    tmp=strsplit(output,'\n');
    dbPosList{i}=[i*ones(length(tmp)-1,1) str2double(tmp(1:end-1)')]; 
end
dbPosList=cell2mat(dbPosList');


%%% use bamList and sampleNames for original mat rather than yaml
if exist('bamList')==0
    bamList{1}={'NA','NA'};
end
inputParam.sampleNames=strjoin(sampleNames,',');
inputParam.bamPaths=strjoin(bamList{1},',');
message=['imported data at: ' char(datetime('now'))]

p=profile('info');
%writetable(p.FunctionTable,[inputParam.outName 'prof.txt']);
profile resume -history
save([inputParam.outMat],'-v7.3');
message=['saved data at: ' char(datetime('now'))]

%pool=gcp('nocreate');
%if isempty(pool)
%    delete(gcp('nocreate'));
%    pc = parcluster('local');
%    pc.NumWorkers = inputParam.numCPU;
%    parpool(pc, pc.NumWorkers);
%end
%%%find most frequent A and B allele accross all samples, and get read
%%%counts for those alleles
% countsAll=getCounts(Tcell, inputParam);
% for j=1:length(Tcell)
%     [lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[countsAll.Chr countsAll.Pos],'rows');
%     Tcell{j}.RefComb=countsAll.Ref(locb);
%     Tcell{j}.Acomb=countsAll.A(locb);
%     Tcell{j}.Bcomb=countsAll.B(locb);
%     Tcell{j}.AcountsComb=countsAll.Acounts(locb,j);
%     Tcell{j}.BcountsComb=countsAll.Bcounts(locb,j);
%     Tcell{j}.ApopAFcomb=countsAll.ApopAF(locb);
%     Tcell{j}.BpopAFcomb=countsAll.BpopAF(locb);
%     Tcell{j}.CosmicCount=countsAll.cosmicCount(locb);
%end
exonCounts=arrayfun(@(x) sum(Ecell{1}.Chr==x),0:length(inputParam.chrList)-1,'UniformOutput',0);
exonOffsets=cumsum([exonCounts{:}])+1;
% for j=1:length(Tcell)
%     %[lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[countsAll.Chr countsAll.Pos],'rows');
%     Tcell{j}.RefComb=Tcell{j}.Ref;
%     Tcell{j}.Acomb=Tcell{j}.A;
%     Tcell{j}.Bcomb=Tcell{j}.B;
%     Tcell{j}.exonIdxAll=Tcell{j}.exonIdx;
%     for i=1:length(exonOffsets)
%         Tcell{j}.exonIdxAll(Tcell{j}.Chr==i)=Tcell{j}.exonIdx(Tcell{j}.Chr==i)+exonOffsets(i);
%     end
%     Tcell{j}.AcountsComb=Tcell{j}.ACountF+Tcell{j}.ACountR;
%     Tcell{j}.BcountsComb=Tcell{j}.BCountF+Tcell{j}.BCountR;
%     Tcell{j}.BcountsComb(Tcell{j}.A==Tcell{j}.B)=0;
%     Tcell{j}.ApopAFcomb=Tcell{j}.ApopAF;
%     Tcell{j}.BpopAFcomb=Tcell{j}.BpopAF;
% end
posTable.exonIdxAll=posTable.BedIndex_0;
for i=1:length(exonOffsets)
    posTable.exonIdxAll(posTable.Chr_0==i)=posTable.BedIndex_0(posTable.Chr_0==i)+exonOffsets(i);
end

%%%Initial Quality Filtering

P=posTable(:,cellstr({"Chr_0","Pos_0"}));
for i=1:inputParam.sampleCount   
    posTable.(strcat('AcountsComb_',num2str(i-1)))=posTable.(strcat('ACountF_',num2str(i-1)))+posTable.(strcat('ACountR_',num2str(i-1)));
    posTable.(strcat('BcountsComb_',num2str(i-1)))=posTable.(strcat('BCountF_',num2str(i-1)))+posTable.(strcat('BCountR_',num2str(i-1)));
    posTable.(strcat('BcountsComb_',num2str(i-1)))(posTable.A_0==posTable.B_0)=0;
end
if inputParam.NormalSample>0
    posTable.aIdx=posTable.(strcat('AcountsComb_',num2str(inputParam.NormalSample)))<posTable.(strcat('BcountsComb_',num2str(inputParam.NormalSample)));
else
    posTable.aIdx=posTable.BpopAF_0>posTable.ApopAF_0;
end
for i=1:inputParam.sampleCount
    T=posTable(:,~cellfun('isempty',regexp(posTable.Properties.VariableNames,strcat('_',num2str(i-1)))) | cellfun('isempty',regexp(posTable.Properties.VariableNames,'_')));    
    T.Properties.VariableNames=regexprep(T.Properties.VariableNames,strcat('_',num2str(i-1)),'');
    E=Ecell{i};
    idx1=T.ApopAF+T.BpopAF>1 & T.ApopAF<T.BpopAF;
    T.BpopAF(idx1)=1-T.ApopAF(idx1)-3*inputParam.pvFreqIndel;
    idx2=T.ApopAF+T.BpopAF>1 & T.ApopAF>=T.BpopAF;
    T.ApopAF(idx2)=1-T.BpopAF(idx2)-3*inputParam.pvFreqIndel;
    Tcell{i}=T;
    homPos=(T.BcountsComb<inputParam.minBCount & ~T.aIdx) | (T.AcountsComb<inputParam.minBCount & T.aIdx);
    [~,postTrust]=qualDiscrimAutoCut(T,E,homPos,inputParam);
    if i==1
        P.trust=postTrust(:,2);
    else
        P.trust=[P.trust  postTrust(:,2)];
    end
    clear T E;
end
filtPos=max(P.trust,[],2)>inputParam.pGoodThresh;
message=['initial quality filtering at: ' char(datetime('now'))]

%%%Find likely het positions
tIdx=setdiff(1:inputParam.sampleCount,inputParam.NormalSample);
if inputParam.NormalSample>0
    commonHet=min(posTable{:,~cellfun('isempty',regexp(posTable.Properties.VariableNames,'popAF'))},[],2)>inputParam.minHetPopFreq & posTable.(strcat('BcountsComb_',num2str(inputParam.NormalSample)))>=inputParam.minBCount;
else
    commonHet=min(posTable{:,~cellfun('isempty',regexp(posTable.Properties.VariableNames,'popAF'))},[],2)>inputParam.minHetPopFreq;
end
posTable.hetPos=commonHet & filtPos;
message=['preliminary variant classification at: ' char(datetime('now'))]
p=profile('info');
%writetable(p.FunctionTable,[inputParam.outName 'prof.txt']);
profile resume -history
save([inputParam.outMat],'-v7.3');
gather(head(posTable))
%%% Filters exon data and segment
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
    hetData{i}=gather(posTable(posTable.hetPos,~cellfun('isempty',regexp(posTable.Properties.VariableNames,strcat('_',num2str(i-1)))) | cellfun('isempty',regexp(posTable.Properties.VariableNames,'_'))));    
    message=['gathered hets at: ' char(datetime('now'))]
    hetData{i}.Properties.VariableNames=regexprep(hetData{i}.Properties.VariableNames,strcat('_',num2str(i-1)),'');
    [segs{i},bafSegs{i}]=segmentData(exonRD{i},hetData{i},inputParam.cnaAlpha);
    clear E;
end
p=profile('info');
%writetable(p.FunctionTable,[inputParam.outName 'prof.txt']);
profile resume -history
save([inputParam.outMat],'-v7.3');
%%%Merge segment boundaries across data types and samples
segsMerged=mergeSegments([segs bafSegs],exonRD,hetData,inputParam);

%%%Make sure segments extend to ends of chromosome
chrNum=unique(Ecell{1}.Chr);
for j=1:length(chrNum)
    i=chrNum(j);
    idx1=find(segsMerged(:,1)==i,1,'first');
    idx2=find(segsMerged(:,1)==i,1,'last');
    if(isempty(idx1))
        segsMerged=[segsMerged; NaN(1,4)];
        segsMerged(end,1)=i;
        idx1=size(segsMerged,1);
        idx2=idx1;
    end
    segsMerged(idx1,2)=min(Ecell{1}{Ecell{1}{:,1}==i,2});
    %segsMerged(idx1,2)=1;
    segsMerged(idx2,3)=max(Ecell{1}{Ecell{1}{:,1}==i,3});
    %segsMerged(idx2,3)=2.5E8;
end
message=['segmented data at: ' char(datetime('now'))]

%%%count database positions in segments
idxDB=getPosInRegions(dbPosList,segsMerged(:,1:3));
dbCounts=hist(idxDB,1:size(segsMerged,1))';

%%%set initial sample fractions
if inputParam.NormalSample<1
    f=[inputParam.priorF(tIdx) diag(inputParam.priorF(tIdx)-0.05)+0.05 (diag(inputParam.priorF(tIdx)-0.05)+0.05)*0.5];
else
    f=inputParam.priorF(tIdx);
end
inputParam.numClones=size(f,2);

%%%Call Copy Number 
%posTable.idxSeg=arrayfun(@(x) getPosInRegions([x.Chr_0 x.Pos_0],segsMerged{:,1:3}),posTable);
posTable.idxSeg=min(posTable.Chr_0,0);
for i=1:size(segsMerged,1)
    posTable.idxSeg(posTable.Chr_0==segsMerged(i,1) & posTable.Pos_0>= segsMerged(i,2) & posTable.Pos_0<=segsMerged(i,3))=i;
end

for i=1:inputParam.sampleCount
    diploidPos=hetData{i}.BcountsComb./hetData{i}.ReadDepthPass>inputParam.minHetAF;
    cInit(i,:)=nanmedian(2*hetData{i}.ControlRD(diploidPos)./hetData{i}.ReadDepthPass(diploidPos));
    wInit(i,:)=nanmedian(exonRD{i}(:,4));
end
dbPos=filtPos & min(posTable{:,~cellfun('isempty',regexp(posTable.Properties.VariableNames,'popAF'))},[],2)>inputParam.maxSomPopFreq;
dbData=gather(posTable(dbPos,:));
message=['gathered dbPos at: ' char(datetime('now'))]
hetData{1}.idxSeg=getPosInRegions([hetData{1}.Chr hetData{1}.Pos],segsMerged(:,1:3));
[N, M, Ftable, log2FC, cnaIdx]=callCNAmulti(hetData,exonRD,segsMerged,inputParam,[cInit(:);wInit(:);f(:)],dbData,dbCounts);
segsTable=array2table(segsMerged(:,1:3),'VariableNames',{'Chr','StartPos','EndPos'});
segsTable.N=N;
segsTable.M=M;
segsTable.F=Ftable;
segsTable.cnaIdx=cnaIdx;
for i=1:length(Tcell)
    segsTable.F(:,i)=Ftable(:,i);
    segsTable.W(:,i)=wInit(i);
    segsTable.log2FC(:,i)=log2FC(:,i);
end
p=profile('info');
%writetable(p.FunctionTable,[inputParam.outName 'prof.txt']);
profile resume -history
save([inputParam.outMat],'-v7.3');

%%%Add copy number info to data table
%Ecell{1}.idxSegStart=getPosInRegionsSplit([Ecell{1}.Chr Ecell{1}.StartPos],segsTable{:,1:3},inputParam.blockSize);
posTable.N=arrayfun(@(x) segsTable.N(x),posTable.idxSeg);
posTable.M=arrayfun(@(x) segsTable.M(x),posTable.idxSeg);
for i=1:inputParam.sampleCount
    posTable.(strcat('cnaF_',num2str(i-1)))=arrayfun(@(x) segsTable.F(x,i),posTable.idxSeg);
    posTable.(strcat('BmeanBQ_',num2str(i-1)))(posTable.(strcat('BcountsComb_',num2str(i-1)))==0)=inputParam.defaultBQ;
    posTable.(strcat('W_',num2str(i-1)))=arrayfun(@(x) segsTable.W(x,i),posTable.idxSeg);
end
message=['initial copy number calls: ' char(datetime('now'))]
gather(head(posTable));
%%%repeat fitting and variant calling until converges
i=1;
somPosOld=min(posTable.Chr_0,0);
W=wInit;
CNAscale=cInit;
indelPos=(posTable.A_0>4 | posTable.B_0>4);
while(true)
    %%% if there is a normal sample, find variant type probabilities
    %%% compared to the normal
    if inputParam.NormalSample>0
        n=1;
        NormalSample=inputParam.NormalSample;
        inputParam.NormalSample=1;
        inputParam.numClones=2;
        sampleCount=inputParam.sampleCount;
        inputParam.sampleCount=2;
        for j=1:inputParam.sampleCount
            if j~=NormalSample
                idx=~cellfun('isempty',regexp(posTable.Properties.VariableNames,strcat('_',num2str(i-1)))) | ~cellfun('isempty',regexp(posTable.Properties.VariableNames,strcat('_',num2str(inputParam.NormalSample)))) | cellfun('isempty',regexp(posTable.Properties.VariableNames,'_'))
                T=posTable(:,idx);
                T.Properties.VariableNames=regexprep(T.Properties.VariableNames,strcat('_',num2str(j-1)),'_1');
                [postCombPair{n},~,pDataCombPair{n},~,~,alleleIdPair{n}]=jointSNV(T, [max(f(n,:)) max(f(n,:))./2],W([inputParam.NormalSample j]), inputParam);
                pSomaticPair=postCombPair{n}.Somatic;
                n=n+1;
            else
                pSomaticPair=min(posTable.Chr_0,0);
            end
            if j==1
                P.SomaticPair=pSomaticPair;
            else
                P.SomaticPair=[P.SomaticPair pSomaticPair];
            end
        end
        inputParam.numClones=size(f,2);
        inputParam.NormalSample=NormalSample;
        inputParam.sampleCount=sampleCount;
    else
        P.SomaticPair=min(P.trust,0);
    end
    %%% jointly find probabilities of variant type
    [postComb, ~,clones,prior,alleleId]=jointSNV(posTable, f, W, inputParam);
    %%% add call probabilities to table and find quality scores
    if i==1
        P=[P postComb];
    else
        P.Somatic=postComb.Somatic;
        P.Het=postComb.Het;
        P.Hom=postComb.Hom;
        P.NonDip=postComb.NonDip;
        P.priorSomatic=prior.Somatic;
        P.priorHet=prior.Het;
        P.priorHom=prior.Hom;
        P.priorNonDip=prior.nonDiploid;
    end
    for j=1:inputParam.sampleCount
        T=posTable(:,~cellfun('isempty',regexp(posTable.Properties.VariableNames,strcat('_',num2str(j-1)))) | cellfun('isempty',regexp(posTable.Properties.VariableNames,'_')));
        T.Properties.VariableNames=regexprep(T.Properties.VariableNames,strcat('_',num2str(j-1)),'');
        E=Ecell{j};
        homPos=P.Hom>0.5  & P.SomaticPair(:,j)<0.5 | (P.Somatic>0.5 & inputParam.NormalSample==j) | (T.BcountsComb<inputParam.minBCount & ~T.aIdx) | (T.AcountsComb<inputParam.minBCount & T.aIdx);
        if j==1
            P.homPos=homPos;
        else
            P.homPos=[P.homPos homPos];
        end
        [F{j},postTrust]=qualDiscrimAutoCut(T,E,homPos,inputParam);
        P.trust(:,j)=postTrust(:,2);
        %P.artifact(:,j)=postArtifact(:,1);
        clear T E;
    end
    %%% classify variants based on type probabilites and quality scores
    [Filter,filtPos]=callVariants(posTable,P,inputParam);
    hetPos=Filter=='GermlineHetPASS';
    somPos=Filter=='SomaticPASS' | Filter=='SomaticPairPASS';
    message=['called variants iteration: ' num2str(i) ' at ' char(datetime('now'))]
    ['Total Somatic positions: ' num2str(sum(somPos))]
    ['SomaticPair positions: ' num2str(sum(strcmp(Filter,'SomaticPairPASS')))]
    ['Het positions: ' num2str(sum(hetPos))]
    p=profile('info');
    %writetable(p.FunctionTable,[inputParam.outName 'prof.txt']);
    profile resume -history
    save([inputParam.outMat],'-v7.3');
    %%% check for convergence
    if sum(abs(somPos-somPosOld))/sum(somPos)<=0.05 || i>=inputParam.maxIter
        break;
    end
    for j=1:inputParam.sampleCount
           hetData{j}=gather(posTable(hetPos,~cellfun('isempty',regexp(posTable.Properties.VariableNames,strcat('_',num2str(j-1)))) | cellfun('isempty',regexp(posTable.Properties.VariableNames,'_'))));    
           hetData{j}.Properties.VariableNames=regexprep(hetData{j}.Properties.VariableNames,strcat('_',num2str(j-1)),'');
    end
    dbPos=filtPos & min(posTable{:,~cellfun('isempty',regexp(posTable.Properties.VariableNames,'popAF'))},[],2)>inputParam.maxSomPopFreq;
    dbData=gather(posTable(dbPos,:));
    somData=gather(posTable(somPos,:));
    message=['gathered positions at: ' char(datetime('now'))]
    [segsTable, W, f, CNAscale, nll, t{i}]=fitCNAmulti(hetData,somData,dbData,exonRD,segsMerged,inputParam,f,CNAscale,W,dbCounts)
    for j=1:length(Tcell)
        T=Tcell{j};
        idx=getPosInRegionSplit([T.Chr T.Pos],segsTable{:,1:3},inputParam.blockSize);
        T.NumCopies=segsTable.N(idx);
        T.MinAlCopies=segsTable.M(idx);
        T.cnaF=segsTable.F(idx,j);
        Tcell{j}=T;
    end
    message=['fit cna iteration: ' num2str(i) ' at ' char(datetime('now'))]
    p=profile('info');
    %writetable(p.FunctionTable,[inputParam.outName 'prof.txt']);
    profile resume -history
    save([inputParam.outMat],'-v7.3');
    i=i+1
    somPosOld=somPos;
    inputParam.numClones=size(f,2);   
end
['converged at iteration' num2str(i)]

%%% sort sample fractions and reassign clone numbers
[~,ord]=sort(mean(f,1),2,'descend');
fSort=f(:,ord);
[~,cloneIdSort]=ismember(clones,ord);
ord=[ord size(f,2)+1];
[~,segsTable.cnaIdx]=ismember(segsTable.cnaIdx,ord);
p=profile('info');
%writetable(p.FunctionTable,[inputParam.outName 'prof.txt']);
profile resume -history
save([inputParam.outMat],'-v7.3');
%%% classify variants based on type probabilities and quality scores
[Filter,~,somaticDetected,trustScore]=callVariants(posTable,P,inputParam);

%%%Get correct alleleIds for "SomaticPair" variants
if inputParam.NormalSample>0
    [~,mIdx]=max(P.SomaticPair,[],2);
    pairPos=strncmp(Filter,'SomaticPair',11);
    if sum(pairPos)>0
        for i=1:length(alleleIdPair)
            alleleId(pairPos & mIdx==i)=alleleIdPair{i}(pairPos & mIdx==i);
        end
    end
end


%%%write output files
filterList=categories(Filter);
somCalls=filterList(startsWith(filterList,'Somatic'));
somPosAll=ismember(Filter,somCalls);

for j=1:inputParam.sampleCount
    somDataAll{j}=gather(posTable(somPosAll,~cellfun('isempty',regexp(posTable.Properties.VariableNames,strcat('_',num2str(j-1)))) | cellfun('isempty',regexp(posTable.Properties.VariableNames,'_'))));
    somDataAll{j}.Properties.VariableNames=regexprep(somDataAll{j}.Properties.VariableNames,strcat('_',num2str(j-1)),'');
end


sampleFrac=writeJointVCF(somDataAll,gather(P(somPosAll,:)),fSort,cloneIdSort(somPosAll,:),alleleId(somPosAll,:),gather(cellstr(Filter(somPosAll,:))),somaticDetected(somPosAll,:),gather(trustScore(somPosAll,:)),inputParam,'somatic.all');
segsTableCond=writeSegVCF(segsTable,exonRD,CNAscale,hetData,inputParam);
p=profile('info');
%writetable(p.FunctionTable,[inputParam.outName 'prof.txt']);
profile resume -history
save([inputParam.outMat],'-v7.3');
writeCloneSummary(segsTable,exonRD,somDataAll,fSort,gather(cloneIdSort(somPosAll,:)),inputParam,gather(cellstr(Filter(somPosAll,:))),gather(somaticDetected(somPosAll,:)),sampleFrac);

plotPos=hetPos | somPos;
for j=1:inputParam.sampleCount
    plotData{j}=gather(posTable(plotPos,~cellfun('isempty',regexp(posTable.Properties.VariableNames,strcat('_',num2str(j-1)))) | cellfun('isempty',regexp(posTable.Properties.VariableNames,'_'))));
    plotData{j}.Properties.VariableNames=regexprep(plotData{j}.Properties.VariableNames,strcat('_',num2str(j-1)),'');
end
plotCNAandVAF(exonRD,segsTable,CNAscale,fSort,plotData,somPos(plotPos),hetPos(plotPos),cloneIdSort(plotPos),inputParam);

% for i=1:length(F)
%     writetable([Tcell{i}(:,1:2) F{i}],[inputParam.outName '_' sampleNames{i} '.qualMetrics.tsv'],'Delimiter','\t','FileType','text','WriteVariableNames',1);
% end

close all;
output=0;
exit;
