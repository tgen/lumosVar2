function [NsegMax, MsegMax, Fout, log2FC, cnaIdx, nll] = callCNAmulti(hetData,exonRD,segsMerged,inputParam,param,dbData,dbCounts)
%callCNA - determine most likley copy number state for each segment
%
% Syntax: [NsegMax, MsegMax, Fout, log2FC, cnaIdx, nll] = callCNAmulti(hetPos,Tcell,exonRD,segsMerged,inputParam,param,dbPos,dbCounts)
%
% Inputs:
%   dataHet: data for germline heterozygous positions with columns:
%       1-'Chr',2-'Pos',3-'ControlRD',4-'TumorRD',5-'Bcount'
%   exonRD: cell array of matrices with one matrix per sample with 
%       the following columns of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segs: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%   inputParam: structure with fields: cnaPrior, minLik
%   param: vector of paramaters of length 2*numClones+1
%       param(1:numSamples) - copy number scaling constant
%       param(numSamples+1:2*numSamples)=W (controls width of allele frequency dist)
%       param(2*numSamples+1:end)=f (sample fractions)
%   dbPos: logical vector indicating positions in the Tcell tables that are
%       variants in population databases
%   dbCounts: counts of total number of possible database positions in each
%       segment
%   
% Outputs:
%   NsegMax: total copies by segment
%   MsegMax: minor allele copies by segment
%   Fout: sampleFraction by segment
%   log2FC: log2 fold change of tumor/control
%   cnaIdx: index of clonal group for each segment
%   nll: sum of segment likliehoods
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain, fitCNAmulti

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 7-May-2018

%------------- BEGIN CODE --------------

%%% read inputs
D=table();
D.Chr=hetData{1}.Chr;
D.Pos=hetData{1}.Pos;
E=table();
E.Chr=exonRD{1}(:,1);
E.StartPos=exonRD{1}(:,2);
E.EndPos=exonRD{1}(:,3);
for i=1:inputParam.sampleCount
    D.ExpReadCount(:,i)=hetData{i}.ControlRD;
    D.TotalReadCount(:,i)=max(hetData{i}.ReadDepthPass,hetData{i}.AcountsComb+hetData{i}.BcountsComb);
    D.MinorReadCount(:,i)=hetData{i}.BcountsComb;
    E.TumorRD(:,i)=exonRD{i}(:,4);
    E.NormalRD(:,i)=exonRD{i}(:,5);
end
CNAscale=param(1:inputParam.sampleCount);
W=param(inputParam.sampleCount+1:2*inputParam.sampleCount);
tIdx=setdiff(1:inputParam.sampleCount,inputParam.NormalSample);
if inputParam.NormalSample>0
    f=[zeros(inputParam.sampleCount,inputParam.numClones) ones(inputParam.sampleCount,1)];
    f(tIdx,1:end-1)=reshape(param(2*inputParam.sampleCount+1:end),[],inputParam.numClones);
else
    f=[reshape(param(2*inputParam.sampleCount+1:end),[],inputParam.numClones) ones(inputParam.sampleCount,1)];
end
if ~isfield(inputParam,'contamIdx')
    inputParam.contamIdx=[];
end

%%% find means accross segments
meanTumorRDexon=NaN(size(segsMerged,1),inputParam.sampleCount);
meanNormalRDexon=NaN(size(segsMerged,1),inputParam.sampleCount);
meanTumorRD=NaN(size(segsMerged,1),inputParam.sampleCount);
meanMinorRD=NaN(size(segsMerged,1),inputParam.sampleCount);
for i=1:inputParam.sampleCount
    meanTumorRDexon(:,i)=getMeanInRegions([E.Chr E.StartPos],E.TumorRD(:,i),segsMerged)+1;
    meanNormalRDexon(:,i)=getMeanInRegions([E.Chr E.StartPos],E.NormalRD(:,i),segsMerged)+1;
    meanTumorRD(:,i)=getMeanInRegions([D.Chr D.Pos],D.TotalReadCount(:,i),segsMerged);
    if ismember(i,inputParam.contamIdx)
        meanMinorRD(:,i)=nan(size(segsMerged,1),1);
    else
        meanMinorRD(:,i)=getMeanInRegions([D.Chr D.Pos],D.MinorReadCount(:,i),segsMerged);
    end
end
meanTumorRDexon(isnan(meanTumorRDexon))=0;

%%% calculate log2FC
log2FC=NaN(size(segsMerged,1),inputParam.sampleCount);
for i=1:inputParam.sampleCount
    log2FC(:,i)=log2((CNAscale(i)./2).*meanTumorRDexon(:,i)./(meanNormalRDexon(:,i)));
end

%%% find copy number per segment and clone
NsegSample=NaN(size(segsMerged,1),size(f,2),inputParam.sampleCount);
for i=1:size(f,2)
    for j=1:inputParam.sampleCount
        if(f(j,i)==0)
            NsegSample(:,i,j)=2*ones(size(segsMerged,1),1);
        else
            NsegSample(:,i,j)=max((CNAscale(j)*(meanTumorRDexon(:,j)./(meanNormalRDexon(:,j)))-2*(1-f(j,i)))/f(j,i),0);
         end
    end
end
NsegSample(NsegSample<0)=0;
Nseg=NaN(size(segsMerged,1),size(f,2),2);
for i=1:size(f,2)
    Nseg(:,i,1)=floor(squeeze(NsegSample(:,i,:))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
    Nseg(:,i,2)=Nseg(:,i,1)+1;
end

%%% find minor allele copy number per segment and clone
MsegSample=NaN(size(segsMerged,1),size(f,2),inputParam.sampleCount,2);
for i=1:size(f,2)
    for j=1:inputParam.sampleCount
        if(f(j,i)==0)
            MsegSample(:,i,j,1)=ones(size(segsMerged,1),1);
            MsegSample(:,i,j,2)=ones(size(segsMerged,1),1);
        else
            fIdx=isfinite(meanMinorRD(:,j));
            MsegSample(fIdx,i,j,1)=Nseg(fIdx,i,1).*meanMinorRD(fIdx,j)./meanTumorRD(fIdx,j)+2.*(meanMinorRD(fIdx,j)./meanTumorRD(fIdx,j)).*(1-f(j,i))./f(j,i)-(1-f(j,i))./f(j,i);
            MsegSample(fIdx,i,j,2)=Nseg(fIdx,i,2).*meanMinorRD(fIdx,j)./meanTumorRD(fIdx,j)+2.*(meanMinorRD(fIdx,j)./meanTumorRD(fIdx,j)).*(1-f(j,i))./f(j,i)-(1-f(j,i))./f(j,i);
            MsegSample(~isfinite(meanMinorRD(:,j)),i,j,1)=min(Nseg(~isfinite(meanMinorRD(:,j)),i,1)-1,1);
            MsegSample(~isfinite(meanMinorRD(:,j)),i,j,2)=min(Nseg(~isfinite(meanMinorRD(:,j)),i,2)-1,1);
        end
    end
    MsegSample(:,i,j,1)=min(MsegSample(:,i,j,1),Nseg(:,i,1)-1);
    MsegSample(:,i,j,2)=min(MsegSample(:,i,j,2),Nseg(:,i,2)-1);
end
MsegSample(MsegSample<0)=0;
Mseg=NaN(size(segsMerged,1),size(f,2),2,2);
for i=1:size(f,2)
    Mseg(:,i,1,1)=min(floor(squeeze(MsegSample(:,i,:,1))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i))),max(Nseg(:,i,1)-1,0));
    Mseg(:,i,1,2)=min(ceil(squeeze(MsegSample(:,i,:,1))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i))),max(Nseg(:,i,1)-1,0));
    Mseg(:,i,2,1)=min(floor(squeeze(MsegSample(:,i,:,2))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i))),max(Nseg(:,i,2)-1,0));
    Mseg(:,i,2,2)=min(ceil(squeeze(MsegSample(:,i,:,2))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i))),max(Nseg(:,i,2)-1,0));
end

%%%find germline priors
germPrior=inputParam.priorGermCNV*ones(size(Nseg,1),2,2);
sexChr=regexp(inputParam.sexChr,',','split');
if cellfun('length',(regexp('',',','split')))==0
    chrList=cellstr(num2str(inputParam.autosomes','%-d'));
else
    chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr'];
end
if cellfun('length',(regexp('',',','split')))>0
    for i=1:length(sexChr)
        chrIdx=find(strcmp(sexChr(i),chrList));
        lia=ismember(Nseg(:,end,1),inputParam.(sexChr{i}));
        germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,1)-Mseg(:,end,1,1)<=1,1,1)=0.5;
        germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,1)-Mseg(:,end,1,2)<=1,1,2)=0.5;
        lia=ismember(Nseg(:,end,2),inputParam.(sexChr{i}));
        germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,2)-Mseg(:,end,2,1)<=1,2,1)=0.5;
        germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,2)-Mseg(:,end,2,2)<=1,2,2)=0.5;
    end
end
    

%%% lookup copy number for positions and exons
idx=hetData{1}.idxSeg;
Nmat=Nseg(idx,:,:);
Mmat=Mseg(idx,:,:,:);
germPriorMat=germPrior(idx,:,:);
idxExon=getPosInRegions([E.Chr E.StartPos],segsMerged);
NmatExon=Nseg(idxExon,:,:);
germPriorExon=germPrior(idxExon,:,:);
%idxAll=getPosInRegions([Tcell{1}.Chr Tcell{1}.Pos], segsMerged);


%%% find prior of copy number
priorCNA=nan(size(NmatExon));
for i=1:length(inputParam.cnaPrior)-1
    priorCNA(NmatExon==i-1)=inputParam.cnaPrior(i);
end
priorCNA(NmatExon>=length(inputParam.cnaPrior)-1)=inputParam.cnaPrior(end);
priorCNA=repmat(priorCNA,1,1,1,2);
priorCNA(:,end,:,:)=germPriorExon;
priorMinAllele=nan(size(Mmat));
for i=1:length(inputParam.minAllelePrior)-1
    priorMinAllele(Mmat==i-1)=inputParam.minAllelePrior(i);
end
priorMinAllele(Mmat>=length(inputParam.minAllelePrior)-1)=inputParam.minAllelePrior(end);
priorMinAllele(:,end,:,:)=germPriorMat;

%%% find prior of clonal group sample fractions
posCounts=accumarray(idxExon,E.EndPos-E.StartPos,[size(segsMerged,1) 1],[],NaN);
pHet=height(hetData{1})./sum(posCounts);
priorCNAf=NaN(size(Mmat));
fDiff=nan(size(f,2)-1,1);
for i=1:size(f,2)-1
    if size(f,1)>1
        fDiff(i)=geomean(pdist(f(:,i),'cityblock')+0.05);
    else
        fDiff(i)=geomean(pdist([0; f(:,i); 1],'cityblock')+0.05);
    end
end
if length(inputParam.priorF)>1
    priorFDiff=geomean(pdist(inputParam.priorF,'cityblock')+0.05);
else
    priorFDiff=geomean(pdist([0; inputParam.priorF; 1],'cityblock')+0.05);
end

%%%find number of somatic and germline variants called at database
%%%positions
corrSeg=NaN(size(Nseg,1),size(f,2),inputParam.sampleCount,2,2);
hetCountOrig=hist(idx,1:size(segsMerged,1))';
%hetPosDB=zeros(length(hetPos),1,size(f,2),2);
exonCounts=hist(idxExon,1:size(segsMerged,1))';
idxDB=dbData.idxSeg;
%idxDB=getPosInRegionSplit([Tcell{1}.Chr(dbPos) Tcell{1}.Pos(dbPos)],segsMerged(:,1:3),inputParam.blockSize);
hetCount=NaN(size(Mseg));
somDBcounts=NaN(size(Mseg));
dbLik=NaN(size(Mseg));
%Tdb=cell(size(Tcell));
if inputParam.NormalSample<1
    for i=1:size(f,2)
        for k=1:2
            for m=1:2
                dbData.N=Nseg(idxDB,i,k);
                dbData.M=Mseg(idxDB,i,k);
                for j=1:inputParam.sampleCount
                    dbData.(strcat('cnaF_',num2str(j-1)))=f(j,i)*ones(height(dbData),1);
		    dbData.(strcat('W_',num2str(j-1)))=W(j)*ones(height(dbData),1);
                end
                postComb=jointSNV(dbData, f(tIdx,1:end-1), W, inputParam);
                somDBpos=postComb.Somatic>inputParam.pSomaticThresh;
                hetPosDB=postComb.Het>inputParam.pGermlineThresh;
                %[lia,locb]=ismember([Tcell{1}.Chr Tcell{1}.Pos],[postComb.Chr postComb.Pos],'rows');
                %hetPosDB(lia,i,k,m)=postComb.Het(locb(lia))>inputParam.pGermlineThresh;
                hetCount(:,i,k,m)=hist(idxDB(hetPosDB==1),1:size(segsMerged),1);
                somDBcounts(:,i,k,m)=hist(idxDB(somDBpos),1:size(segsMerged,1));
                dbLik(:,i,k,m)=1-binocdf(somDBcounts(:,i,k,m),dbCounts,inputParam.priorSomaticSNV);
            end
        end
    end
end

%%% find likelihoods of read depths and allele fractions
expReadCount=NaN(height(E),size(f,2),inputParam.sampleCount,2);
depthlik=NaN(height(E),size(f,2),inputParam.sampleCount,2);
corr=NaN(size(Mmat));
hetlik=NaN(height(D),size(f,2),inputParam.sampleCount,2,2);
pHetDetect=NaN(size(segsMerged,1),size(f,2),inputParam.sampleCount,2,2);
hetCountLik=NaN(size(segsMerged,1),size(f,2),inputParam.sampleCount,2,2);
hetExp=NaN(size(segsMerged,1),size(f,2),inputParam.sampleCount,2,2);
segLik=NaN(size(segsMerged,1),size(f,2),inputParam.sampleCount,2,2);
for i=1:size(f,2)
    for k=1:2
        for j=1:inputParam.sampleCount
            expReadCount(:,i,j,k)=f(j,i)*E.NormalRD(:,j).*NmatExon(:,i,k)./CNAscale(j)+(1-f(j,i))*E.NormalRD(:,j)*2./CNAscale(j);
            depthlik(:,i,j,k)=poisspdf(round(E.TumorRD(:,j)),round(expReadCount(:,i,j,k)))+inputParam.minLik;
            for m=1:2
                if i==size(f,2)
                    priorCNAf(:,i,k,m)=germPriorMat(:,k,m);
                else
                    priorCNAf(:,i,k,m)=betacdf(fDiff(i),inputParam.alphaF,(inputParam.alphaF-1)./priorFDiff-inputParam.alphaF+2)+inputParam.minLik;
                    priorCNAf(Nmat(:,i,k)==2 & Mmat(:,i,k,m)==1,i,k,m)=1;
                end
                corr(:,i,j,k,m)=(f(j,i).*Mmat(:,i,k,m)+(1-f(j,i)))./(f(j,i).*Nmat(:,i,k)+(1-f(j,i)).*2);
                if f(j,i)==1
                    corr(Nmat(:,i,k)==0,i,j,k,m)=0;
                else
                    corr(Nmat(:,i,k)==0,i,j,k,m)=0.5;
                end
                corr(corr(:,i,j,k,m)<0,i,j,k,m)=0;
                corr(corr(:,i,j,k,m)>1,i,j,k,m)=1;
                corrSeg(:,i,j,k,m)=(f(j,i).*Mseg(:,i,k,m)+(1-f(j,i)))./(f(j,i).*Nseg(:,i,k)+(1-f(j,i)).*2);
                if f(j,i)==1
                    corrSeg(Nseg(:,i,k)==0,i,j,k,m)=0;
                else
                    corrSeg(Nseg(:,i,k)==0,i,j,k,m)=0.5;
                end
                corrSeg(corrSeg(:,i,j,k,m)<0,i,j,k,m)=0;
                corrSeg(corrSeg(:,i,j,k,m)>1,i,j,k,m)=1;
                hetlik(:,i,j,k,m)=bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j)*corr(:,i,j,k,m),W(j)*(1-corr(:,i,j,k,m)))+(corr(:,i,j,k,m)./0.5).*bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j)*(1-corr(:,i,j,k,m)),W(j)*corr(:,i,j,k,m))+inputParam.minLik;
                hetlik(corr(:,i,j,k,m)==0,i,j,k,m)=inputParam.minLik;
                pHetDetect(:,i,j,k,m)=binocdf(inputParam.minBCount,round(meanTumorRDexon(:,j)),corrSeg(:,i,j,k,m),'upper');
                hetCountLik(:,i,j,k,m)=binocdf(hetCount(:,i,k,m),dbCounts,pHetDetect(:,i,j,k,m)*pHet);
                hetExp(:,i,j,k,m)=dbCounts.*pHetDetect(:,i,j,k,m)*pHet;
                if ismember(j,inputParam.contamIdx)
                    hetlik(:,i,j,k,m)=NaN;
                end
                segLik(:,i,j,k,m)=nansum([(hetCountOrig./sum(hetCountOrig)).*getMeanInRegions([D.Chr D.Pos],log(hetlik(:,i,j,k,m))+log(priorCNAf(:,i,k,m))+log(priorMinAllele(:,i,k,m)),segsMerged) (exonCounts./sum(exonCounts)).*getMeanInRegions([E.Chr E.StartPos],log(depthlik(:,i,j,k))+log(priorCNA(:,i,k,m)),segsMerged) (dbCounts./sum(dbCounts)).*log(hetCountLik(:,i,j,k,m))+(posCounts./sum(posCounts)).*log(dbLik(:,i,k,m)+realmin)],2);   
            end
        end
    end
end

%%% find which clone contains most likely CNV per segment
[kMax,kIdx]=max(sum(segLik,3),[],4);
[mMax,mIdx]=max(kMax,[],5);
[segLikMax,cnaIdx]=max(mMax,[],2);
macnIdx=NaN(size(cnaIdx));
cnIdx=NaN(size(cnaIdx));
for i=1:size(f,2)
    macnIdx(cnaIdx==i,:)=mIdx(cnaIdx==i,i);
end
for i=1:size(f,2)
    for m=1:2
        cnIdx(cnaIdx==i & macnIdx==m,:)=squeeze(kIdx(cnaIdx==i & macnIdx==m,i,1,1,m));
    end
end

NsegMax=NaN(size(segsMerged,1),1);
MsegMax=NaN(size(segsMerged,1),1);
for i=1:size(f,2)
    for k=1:2
        NsegMax(cnaIdx==i & cnIdx==k,:)=Nseg(cnaIdx==i & cnIdx==k,i,k);
        for m=1:2
            MsegMax(cnaIdx==i & cnIdx==k & macnIdx==m,:)=Mseg(cnaIdx==i & cnIdx==k & macnIdx==m,i,k,m);
        end
    end
end

Fout=NaN(length(cnaIdx),inputParam.sampleCount);
for j=1:inputParam.sampleCount
    Fout(:,j)=f(j,cnaIdx);
    %Wout(:,j)=W(j,cnaIdx);
end

nll=sum(segLikMax);

return
