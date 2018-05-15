function output=LumosVarMain(paramFile,numCPU,varargin)
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

inputParam=readInputs(paramFile);
cd(inputParam.workingDirectory);

Tcell=cell('');
Ecell=cell('');
%%%% If outmat exits, uses outmat, otherwise generates data tables
if(exist([inputParam.outMat],'file'))
    vars={'T','E','Tcell','Ecell','sampleNames','bamList','inputParam'};
    load([inputParam.outMat],'-mat',vars{:});
    NormalSample=inputParam.NormalSample;
    priorF=inputParam.priorF;
    if exist('sampleNames')==0
        sampleNames=strsplit(inputParam.sampleNames,',');
    end
    inputParam=readInputs(paramFile);
    inputParam.NormalSample=NormalSample;
    if length(priorF)==length(Tcell)
        inputParam.priorF=priorF;
    end
else
    [T, E]=preprocessTumorOnly(inputParam,paramFile);   
    fid=fopen(inputParam.bamList);
    bamList=textscan(fid,'%s');
    sampleCount=length(bamList{1});
    fclose(fid);
    for i=1:sampleCount
        folders=regexp(bamList{1}{i},'/','split');
        names=regexp(folders{end},'\.','split');
        sampleNames{i}=names{1};
        message=['sample ' num2str(i) ': ' sampleNames{i} ' from ' bamList{1}{i}]
    end
    save([inputParam.outMat],'-v7.3');
end

%%%% Initial parallel pool
inputParam.numCPU=str2num(numCPU);
pool=gcp('nocreate');
if isempty(pool)
    delete(gcp('nocreate'));
    pc = parcluster('local');
    pc.NumWorkers = inputParam.numCPU;
    parpool(pc, pc.NumWorkers);
end

%%%% get chromosomes to examine
sexChr=regexp(inputParam.sexChr,',','split');
chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr'];

%%%% find positions within bed file that have a population allele frequency
%%%% greater than maxSomPopFreq
parfor i=1:length(chrList)
    cmd=strcat({'bcftools view -q '},num2str(inputParam.maxSomPopFreq),{' -R '},inputParam.regionsFile,{' '},inputParam.snpVCFpath,chrList(i),inputParam.snpVCFname,' | grep -v ^# | awk ''{ print $2 "\n" }''')
    [~,output]=system(cmd{1});
    tmp=strsplit(output,'\n');
    dbPosList{i}=[i*ones(length(tmp)-1,1) str2double(tmp(1:end-1)')]; 
end
dbPosList=cell2mat(dbPosList');


%%% rename T and E
if isempty(Tcell)
    Tcell=T;
end
if isempty(Ecell)
    Ecell=E;
end

%%% use bamList and sampleNames for original mat rather than yaml
if exist('bamList')==0
    bamList{1}={'NA','NA'};
end
inputParam.sampleNames=strjoin(sampleNames,',');
inputParam.bamPaths=strjoin(bamList{1},',');
message=['imported data at: ' char(datetime('now'))]

save([inputParam.outMat],'-v7.3');
message=['saved data at: ' char(datetime('now'))]

%%%Initial Quality Filtering
P=table();
P.trust=nan(height(Tcell{1}),length(Tcell));
P.artifact=nan(height(Tcell{1}),length(Tcell));
for i=1:size(Tcell,2)
    T=Tcell{i};
    E=Ecell{i};
    idx1=T.ApopAF+T.BpopAF>1 & T.ApopAF<T.BpopAF;
    T.BpopAF(idx1)=1-T.ApopAF(idx1)-3*inputParam.pvFreqIndel;
    idx2=T.ApopAF+T.BpopAF>1 & T.ApopAF>=T.BpopAF;
    T.ApopAF(idx2)=1-T.BpopAF(idx2)-3*inputParam.pvFreqIndel;
    Tcell{i}=T;
    [~,postTrust,postArtifact]=qualDiscrimCalls(T,E,inputParam);
    P.trust(:,i)=postTrust(:,2);
    P.artifact(:,i)=postArtifact(:,1);
    clear T E;
end
filtPos=max(P.trust,[],2)>inputParam.pGoodThresh & max(P.artifact,[],2)<inputParam.pGoodThresh;
message=['initial quality filtering at: ' char(datetime('now'))]

%%%find most frequent A and B allele accross all samples, and get read
%%%counts for those alleles
countsAll=getCounts(Tcell, inputParam);
for j=1:length(Tcell)
    [lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[countsAll.Chr countsAll.Pos],'rows');
    Tcell{j}.RefComb=countsAll.Ref(locb);
    Tcell{j}.Acomb=countsAll.A(locb);
    Tcell{j}.Bcomb=countsAll.B(locb);
    Tcell{j}.AcountsComb=countsAll.Acounts(locb,j);
    Tcell{j}.BcountsComb=countsAll.Bcounts(locb,j);
    Tcell{j}.ApopAFcomb=countsAll.ApopAF(locb);
    Tcell{j}.BpopAFcomb=countsAll.BpopAF(locb);
    Tcell{j}.CosmicCount=countsAll.cosmicCount(locb);
end

%%%Find likely het positions
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
if inputParam.NormalSample>0
    commonHet=(min([Tcell{inputParam.NormalSample}.ApopAFcomb Tcell{inputParam.NormalSample}.BpopAFcomb],[],2)>inputParam.minHetPopFreq) & Tcell{inputParam.NormalSample}.BCountF+Tcell{inputParam.NormalSample}.BCountR>=inputParam.minBCount;
else
    commonHet=min([Tcell{1}.ApopAFcomb Tcell{1}.BpopAFcomb],[],2)>inputParam.minHetPopFreq;
end    
hetPos=commonHet & filtPos;
message=['preliminary variant classification at: ' char(datetime('now'))]


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
    [segs{i},bafSegs{i}]=segmentData(exonRD{i},Tcell{i}(hetPos,:),inputParam.cnaAlpha);
    clear E;
end

%%%Merge segment boundaries across data types and samples
segsMerged=mergeSegments([segs bafSegs],exonRD,Tcell,hetPos,inputParam);

%%%Make sure segments extend to ends of chromosome
numChr=max(Ecell{1}.Chr);
for i=1:numChr
    idx1=find(segsMerged(:,1)==i,1,'first');
    idx2=find(segsMerged(:,1)==i,1,'last');
    if(isempty(idx1))
        segsMerged=[segsMerged; NaN(1,4)];
        segsMerged(end,1)=i;
        idx1=size(segsMerged,1);
        idx2=idx1;
    end
    segsMerged(idx1,2)=min([Tcell{1}.Pos(Tcell{1}.Chr==i); Ecell{1}{Ecell{1}{:,1}==i,2}]);
    segsMerged(idx2,3)=max([Tcell{1}.Pos(Tcell{1}.Chr==i); Ecell{1}{Ecell{1}{:,1}==i,3}]);
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
for i=1:length(Tcell)
    diploidPos=(Tcell{i}.BCountF+Tcell{i}.BCountR)./Tcell{i}.ReadDepthPass>inputParam.minHetAF;
    cInit(i,:)=median(2*Tcell{i}.ControlRD(diploidPos & hetPos)./Tcell{i}.ReadDepthPass(diploidPos & hetPos));
    wInit(i,:)=nanmedian(exonRD{i}(:,4));
end
dbPos=filtPos & min([Tcell{1}.ApopAFcomb Tcell{1}.BpopAFcomb],[],2)>inputParam.maxSomPopFreq;
[N, M, Ftable, log2FC, cnaIdx]=callCNAmulti(hetPos,Tcell,exonRD,segsMerged,inputParam,[cInit(:);wInit(:);f(:)],dbPos,dbCounts);
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

%%%Add copy number info to data table
for i=1:size(Tcell,2)
    T=Tcell{i};
    idx=getPosInRegionSplit([T.Chr T.Pos],segsTable{:,1:3},inputParam.blockSize);
    T.NumCopies=segsTable.N(idx);
    T.MinAlCopies=segsTable.M(idx);
    T.cnaF=segsTable.F(idx,i);
    T.W=segsTable.W(idx);
    T.BmeanBQ(T.BCountF+T.BCountR==0)=inputParam.defaultBQ;
    Tcell{i}=T;
    clear T;
end
message=['initial copy number calls: ' char(datetime('now'))]
save([inputParam.outMat],'-v7.3');

%%%repeat fitting and variant calling until converges
i=1;
somPosOld=zeros(size(hetPos));
W=wInit;
CNAscale=cInit;
indelPos=(Tcell{1}.Acomb>4 | Tcell{1}.Bcomb>4);
while(true)
    %%% if there is a normal sample, find variant type probabilities
    %%% compared to the normal
    if inputParam.NormalSample>0
        n=1;
        NormalSample=inputParam.NormalSample;
        inputParam.NormalSample=1;
        inputParam.numClones=2;
        for j=1:length(Tcell)
            if j~=NormalSample
                idx=[NormalSample j];
                [postCombPair{n},~,pDataCombPair{n},~,~,alleleIdPair{n}]=jointSNV(Tcell(idx), [max(f(n,:)) max(f(n,:))./2],W(idx), inputParam);
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
    %%% jointly find probabilities of variant type
    [postComb, pDataSum(i+1),pDataComb,clones,prior,alleleId]=jointSNV(Tcell, f, W, inputParam);
    %%% add call probabilities to table and find quality scores
    for j=1:length(Tcell)
        [lia,locb]=ismember([Tcell{j}.Chr Tcell{j}.Pos],[postComb.Chr postComb.Pos],'rows');
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
    %%% classify variants based on type probabilites and quality scores
    [Filter,filtPos]=callVariants(Tcell,P,inputParam);
    hetPos=strcmp(Filter,'GermlineHetPASS');
    somPos=strcmp(Filter,'SomaticPASS') | strcmp(Filter,'SomaticPairPASS');
    message=['called variants iteration: ' num2str(i) ' at ' char(datetime('now'))]
    ['Total Somatic positions: ' num2str(sum(somPos))]
    ['SomaticPair positions: ' num2str(sum(strcmp(Filter,'SomaticPairPASS')))]
    ['Het positions: ' num2str(sum(hetPos))]
    save([inputParam.outMat],'-v7.3');
    %%% check for convergence
    if sum(abs(somPos-somPosOld))/sum(somPos)<=0.05 || i>=inputParam.maxIter
        break;
    end
    [segsTable, W, f, CNAscale, nll, t{i}]=fitCNAmulti_v3(hetPos & ~indelPos,somPos & ~indelPos,filtPos & min([Tcell{1}.ApopAFcomb Tcell{1}.BpopAFcomb],[],2)>inputParam.maxSomPopFreq,Tcell,exonRD,segsMerged,inputParam,f,CNAscale,W,dbCounts)
    for j=1:length(Tcell)
        T=Tcell{j};
        idx=getPosInRegionSplit([T.Chr T.Pos],segsTable{:,1:3},inputParam.blockSize);
        T.NumCopies=segsTable.N(idx);
        T.MinAlCopies=segsTable.M(idx);
        T.cnaF=segsTable.F(idx,j);
        Tcell{j}=T;
    end
    message=['fit cna iteration: ' num2str(i) ' at ' char(datetime('now'))]
    save([inputParam.outMat],'-v7.3');
    i=i+1
    somPosOld=somPos;
    inputParam.numClones=size(f,2);   
end
['converged at iteration' num2str(i)]

%%% sort sample fractions and reassign clone numbers
[~,ord]=sort(mean(f,1),2,'descend');
fSort=f(:,ord);
[~,cloneIdSort]=ismember(cloneId,ord);
ord=[ord size(f,2)+1];
[~,segsTable.cnaIdx]=ismember(segsTable.cnaIdx,ord);

save([inputParam.outMat],'-v7.3');
%%% classify variants based on type probabilities and quality scores
[Filter,~,somaticDetected,trustScore,artifactScore]=callVariants(Tcell,P,inputParam);

%%%Get correct alleleIds for "SomaticPair" variants
[~,mIdx]=max(P.SomaticPair,[],2);
pairPos=strncmp(Filter,'SomaticPair',11);
if sum(pairPos)>0
    for i=1:length(alleleIdPair)
        alleleId(pairPos & mIdx==i)=alleleIdPair{i}(pairPos & mIdx==i);
    end
end

%%%write output files
sampleFrac=writeJointVCF(Tcell,P,fSort,cloneIdSort,alleleId,Filter,somaticDetected,trustScore,artifactScore,inputParam);
segsTableCond=writeSegVCF(segsTable,exonRD,CNAscale,Tcell,hetPos,inputParam);
save([inputParam.outMat],'-v7.3');
writeCloneSummary(segsTable,exonRD,Tcell,fSort,cloneIdSort,inputParam,Filter,somaticDetected,sampleFrac);
plotTumorOnly(exonRD,segsTable,CNAscale,fSort,Tcell,somPos,hetPos,cloneIdSort,inputParam);
for i=1:length(F)
    writetable([Tcell{i}(:,1:2) F{i}],[inputParam.outName '_' sampleNames{i} '.qualMetrics.tsv'],'Delimiter','\t','FileType','text','WriteVariableNames',1);
end

close all;
output=0;
exit;
