function [NsegMax, MsegMax, Fout, log2FC, cnaIdx, nll] = callCNAmulti_v3(hetPos,Tcell,exonRD,segsMerged,inputParam,param,dbPos,dbCounts)
%callCNA - determine most likley copy number state for each segment
%
% Syntax: [N, M, F, Wout, log2FC, pCNA] = callCNA(dataHet,exonRD,segs,inputParam,param)
%
% Inputs:
%   dataHet: data for germline heterozygous positions with columns:
%       1-'Chr',2-'Pos',3-'ControlRD',4-'TumorRD',5-'Bcount'
%   exonRD: matrix of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segs: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%   inputParam: structure with fields: cnaPrior, minLik
%   param: vector of paramaters of length 2*numClones+1
%       param(1) - copy number scaling constant
%       param(2:numClones+1)=W (controls width of allele frequency dist)
%       param(numClones+1:2*numClones+1)=f (sample fractions)
%   
% Outputs:
%   N: total copies by segment
%   M: minor allele copies by segment
%   Wout: W parameter by segment
%   log2FC: log2 fold change of tumor/control
%   pCNA: probability of copy number event by exon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, fitCNA, nllCNA

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

%%% read inputs
D=table();
D.Chr=Tcell{1}.Chr(hetPos);
D.Pos=Tcell{1}.Pos(hetPos);
E=table();
E.Chr=exonRD{1}(:,1);
E.StartPos=exonRD{1}(:,2);
E.EndPos=exonRD{1}(:,3);
for i=1:length(Tcell)
    D.ExpReadCount(:,i)=Tcell{i}.ControlRD(hetPos);
    D.TotalReadCount(:,i)=max(Tcell{i}.ReadDepthPass(hetPos),Tcell{i}.AcountsComb(hetPos)+Tcell{i}.BcountsComb(hetPos));
    D.MinorReadCount(:,i)=Tcell{i}.BcountsComb(hetPos);
    E.TumorRD(:,i)=exonRD{i}(:,4);
    E.NormalRD(:,i)=exonRD{i}(:,5);
end
CNAscale=param(1:length(Tcell));
W=param(length(Tcell)+1:2*length(Tcell));
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
if inputParam.NormalSample>0
    f=[zeros(length(Tcell),inputParam.numClones) ones(length(Tcell),1)];
    f(tIdx,1:end-1)=reshape(param(2*length(Tcell)+1:end),[],inputParam.numClones);
else
    f=[reshape(param(2*length(Tcell)+1:end),[],inputParam.numClones) ones(length(Tcell),1)];
end

if ~isfield(inputParam,'contamIdx')
    inputParam.contamIdx=[];
end
%%% find means accross segments
meanTumorRDexon=NaN(size(segsMerged,1),length(Tcell));
meanNormalRDexon=NaN(size(segsMerged,1),length(Tcell));
meanTumorRD=NaN(size(segsMerged,1),length(Tcell));
meanMinorRD=NaN(size(segsMerged,1),length(Tcell));
for i=1:length(Tcell)
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

log2FC=NaN(size(segsMerged,1),length(Tcell));
%%% calculate log2FC
for i=1:length(Tcell)
    log2FC(:,i)=log2((CNAscale(i)./2).*meanTumorRDexon(:,i)./(meanNormalRDexon(:,i)));
end

NsegSample=NaN(size(segsMerged,1),size(f,2),length(Tcell));
%%% find copy number per segment and clone
for i=1:size(f,2)
    for j=1:length(Tcell)
        if(f(j,i)==0)
            NsegSample(:,i,j)=2*ones(size(segsMerged,1),1);
            %MsegSample(:,i,j)=ones(size(segsMerged,1),1);
        else
            NsegSample(:,i,j)=max((CNAscale(j)*(meanTumorRDexon(:,j)./(meanNormalRDexon(:,j)))-2*(1-f(j,i)))/f(j,i),0);
            %MsegSample(isfinite(meanMinorRD(:,j)),i,j)=(NsegSample(isfinite(meanMinorRD(:,j)),i,j)/f(j,i)).*((meanMinorRD(isfinite(meanMinorRD(:,j)),j)./(meanTumorRD(isfinite(meanMinorRD(:,j)),j)))-0.5*(1-f(j,i)));
            %MsegSample(~isfinite(meanMinorRD(:,j)),i,j)=min(NsegSample(~isfinite(meanMinorRD(:,j)),i,j)-1,1);
        end
    end
end
NsegSample(NsegSample<0)=0;
%MsegSample(MsegSample<0)=0;
Nseg=NaN(size(segsMerged,1),size(f,2),2);
for i=1:size(f,2)
    Nseg(:,i,1)=floor(squeeze(NsegSample(:,i,:))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
    Nseg(:,i,2)=Nseg(:,i,1)+1;
    %Mseg(:,i)=round(squeeze(MsegSample(:,i,:))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
end
MsegSample=NaN(size(segsMerged,1),size(f,2),length(Tcell),2);
for i=1:size(f,2)
    for j=1:length(Tcell)
        if(f(j,i)==0)
            %NsegSample(:,i,j)=2*ones(size(segsMerged,1),1);
            MsegSample(:,i,j,1)=ones(size(segsMerged,1),1);
            MsegSample(:,i,j,2)=ones(size(segsMerged,1),1);
        else
            %NsegSample(:,i,j)=max((CNAscale(j)*(meanTumorRDexon(:,j)./(meanNormalRDexon(:,j)))-2*(1-f(j,i)))/f(j,i),0);
            fIdx=isfinite(meanMinorRD(:,j));
            MsegSample(fIdx,i,j,1)=Nseg(fIdx,i,1).*meanMinorRD(fIdx,j)./meanTumorRD(fIdx,j)+2.*(meanMinorRD(fIdx,j)./meanTumorRD(fIdx,j)).*(1-f(j,i))./f(j,i)-(1-f(j,i))./f(j,i);
            MsegSample(fIdx,i,j,2)=Nseg(fIdx,i,2).*meanMinorRD(fIdx,j)./meanTumorRD(fIdx,j)+2.*(meanMinorRD(fIdx,j)./meanTumorRD(fIdx,j)).*(1-f(j,i))./f(j,i)-(1-f(j,i))./f(j,i);
            %MsegSample(isfinite(meanMinorRD(:,j)),i,j,1)=(Nseg(isfinite(meanMinorRD(:,j)),i,1)/f(j,i)).*((meanMinorRD(isfinite(meanMinorRD(:,j)),j)./(meanTumorRD(isfinite(meanMinorRD(:,j)),j)))-0.5*(1-f(j,i)));
            %MsegSample(isfinite(meanMinorRD(:,j)),i,j,2)=(Nseg(isfinite(meanMinorRD(:,j)),i,2)/f(j,i)).*((meanMinorRD(isfinite(meanMinorRD(:,j)),j)./(meanTumorRD(isfinite(meanMinorRD(:,j)),j)))-0.5*(1-f(j,i)));
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
    %Nseg(:,i,1)=floor(squeeze(NsegSample(:,i,:))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
    %Nseg(:,i,2)=Nseg(:,i,1)+1;
    Mseg(:,i,1,1)=min(floor(squeeze(MsegSample(:,i,:,1))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i))),max(Nseg(:,i,1)-1,0));
    Mseg(:,i,1,2)=min(ceil(squeeze(MsegSample(:,i,:,1))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i))),max(Nseg(:,i,1)-1,0));
    Mseg(:,i,2,1)=min(floor(squeeze(MsegSample(:,i,:,2))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i))),max(Nseg(:,i,2)-1,0));
    Mseg(:,i,2,2)=min(ceil(squeeze(MsegSample(:,i,:,2))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i))),max(Nseg(:,i,2)-1,0));
end


germPrior=inputParam.priorGermCNV*ones(size(Nseg,1),2,2);
%[lia,locb]=ismember(segsMerged(:,1),inputParam.autosomes);
%germPrior(lia & Nseg(:,end,1)==2 & Mseg(:,end,1)==1,1)=NaN;
%germPrior(lia & Nseg(:,end,2)==2 & Mseg(:,end,2)==1,2)=NaN;
sexChr=regexp(inputParam.sexChr,',','split');
chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr'];
for i=1:length(sexChr)
    chrIdx=find(strcmp(sexChr(i),chrList));
    lia=ismember(Nseg(:,end,1),inputParam.(sexChr{i}));
    germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,1)-Mseg(:,end,1,1)<=1,1,1)=0.5;
    germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,1)-Mseg(:,end,1,2)<=1,1,2)=0.5;
    lia=ismember(Nseg(:,end,2),inputParam.(sexChr{i}));
    germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,2)-Mseg(:,end,2,1)<=1,2,1)=0.5;
    germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,2)-Mseg(:,end,2,2)<=1,2,2)=0.5;
end
    

%%% lookup copy number for positions and exons
idx=getPosInRegions([D.Chr D.Pos], segsMerged);
Nmat=Nseg(idx,:,:);
Mmat=Mseg(idx,:,:,:);
germPriorMat=germPrior(idx,:,:);
idxExon=getPosInRegions([E.Chr E.StartPos],segsMerged);
NmatExon=Nseg(idxExon,:,:);
germPriorExon=germPrior(idxExon,:,:);
idxAll=getPosInRegions([Tcell{1}.Chr Tcell{1}.Pos], segsMerged);


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

%%% find likelihoods of read counts and depth

posCounts=accumarray(idxExon,E.EndPos-E.StartPos,[size(segsMerged,1) 1],[],NaN);
%pHet=sum(hetPos & dbPos)./sum(dbCounts);
pHet=sum(hetPos)./sum(posCounts);
priorCNAf=NaN(size(Mmat));
%[fMax,fIdx]=max(f(:,1:end-1),[],1);
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

corrSeg=NaN(size(Nseg,1),size(f,2),length(Tcell),2,2);
hetCountOrig=hist(idx,1:size(segsMerged,1))';
%hetPosDB=repmat(hetPos,1,size(f,2),2);
hetPosDB=zeros(length(hetPos),1,size(f,2),2);
exonCounts=hist(idxExon,1:size(segsMerged,1))';
idxDB=getPosInRegionSplit([Tcell{1}.Chr(dbPos) Tcell{1}.Pos(dbPos)],segsMerged(:,1:3),inputParam.blockSize);
%dbCount=hist(idxDB,1:size(segsMerged,1));
hetCount=NaN(size(Mseg));
somDBcounts=NaN(size(Mseg));
dbLik=NaN(size(Mseg));
Tdb=cell(size(Tcell));
if inputParam.NormalSample<1
    for i=1:size(f,2)
        for k=1:2
            for m=1:2
                for j=1:length(Tcell)
                    T=Tcell{j}(dbPos,:);
                    T.NumCopies=Nseg(idxDB,i,k);
                    T.MinAlCopies=Mseg(idxDB,i,k,m);
                    T.cnaF=f(j,i)*ones(height(T),1);
                    T.W=W(j)*ones(height(T),1);
                    Tdb{j}=T;
                end
                postComb=jointSNV_v2(Tdb, f(tIdx,1:end-1), W, inputParam);
                somDBpos=postComb.Somatic>inputParam.pSomaticThresh;
                [lia,locb]=ismember([Tcell{1}.Chr Tcell{1}.Pos],[postComb.Chr postComb.Pos],'rows');
                hetPosDB(lia,i,k,m)=postComb.Het(locb(lia))>inputParam.pGermlineThresh;
                hetCount(:,i,k,m)=hist(idxAll(hetPosDB(:,i,k,m)==1),1:size(segsMerged),1);
                somDBcounts(:,i,k,m)=hist(idxDB(somDBpos),1:size(segsMerged,1));
                dbLik(:,i,k,m)=1-binocdf(somDBcounts(:,i,k,m),dbCounts,inputParam.priorSomaticSNV);
            end
        end
    end
end
expReadCount=NaN(height(E),size(f,2),length(Tcell),2);
depthlik=NaN(height(E),size(f,2),length(Tcell),2);
corr=NaN(size(Mmat));
hetlik=NaN(height(D),size(f,2),length(Tcell),2,2);
pHetDetect=NaN(size(segsMerged,1),size(f,2),length(Tcell),2,2);
hetCountLik=NaN(size(segsMerged,1),size(f,2),length(Tcell),2,2);
hetExp=NaN(size(segsMerged,1),size(f,2),length(Tcell),2,2);
segLik=NaN(size(segsMerged,1),size(f,2),length(Tcell),2,2);
for i=1:size(f,2)
    for k=1:2
        for j=1:length(Tcell)
            expReadCount(:,i,j,k)=f(j,i)*E.NormalRD(:,j).*NmatExon(:,i,k)./CNAscale(j)+(1-f(j,i))*E.NormalRD(:,j)*2./CNAscale(j);
            depthlik(:,i,j,k)=poisspdf(round(E.TumorRD(:,j)),round(expReadCount(:,i,j,k)))+inputParam.minLik;
            for m=1:2
                if i==size(f,2)
                    priorCNAf(:,i,k,m)=germPriorMat(:,k,m);
                else
                    %priorCNAf(:,i,k,m)=betapdf(fMax(i),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF(fIdx(i))-inputParam.alphaF+2)+inputParam.minLik;
                    %priorCNAf(Nmat(:,i,k)==2 & Mmat(:,i,k,m)==1,i,k,m)=betapdf(inputParam.priorF(j),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF(j)-inputParam.alphaF+2);
                    priorCNAf(:,i,k,m)=betacdf(fDiff(i),inputParam.alphaF,(inputParam.alphaF-1)./priorFDiff-inputParam.alphaF+2)+inputParam.minLik;
                    priorCNAf(Nmat(:,i,k)==2 & Mmat(:,i,k,m)==1,i,k,m)=1;
                end
                %corr(:,i,j,k)=f(j,i).*Mmat(:,i,k)./Nmat(:,i,k)+(1-f(j,i))*0.5;
                corr(:,i,j,k,m)=(f(j,i).*Mmat(:,i,k,m)+(1-f(j,i)))./(f(j,i).*Nmat(:,i,k)+(1-f(j,i)).*2);
                if f(j,i)==1
                    corr(Nmat(:,i,k)==0,i,j,k,m)=0;
                else
                    corr(Nmat(:,i,k)==0,i,j,k,m)=0.5;
                end
                corr(corr(:,i,j,k,m)<0,i,j,k,m)=0;
                corr(corr(:,i,j,k,m)>1,i,j,k,m)=1;
                %corrSeg(:,i,j,k)=f(j,i).*Mseg(:,i,k)./Nseg(:,i,k)+(1-f(j,i))*0.5;
                corrSeg(:,i,j,k,m)=(f(j,i).*Mseg(:,i,k,m)+(1-f(j,i)))./(f(j,i).*Nseg(:,i,k)+(1-f(j,i)).*2);
                if f(j,i)==1
                    corrSeg(Nseg(:,i,k)==0,i,j,k,m)=0;
                else
                    corrSeg(Nseg(:,i,k)==0,i,j,k,m)=0.5;
                end
                corrSeg(corrSeg(:,i,j,k,m)<0,i,j,k,m)=0;
                corrSeg(corrSeg(:,i,j,k,m)>1,i,j,k,m)=1;
                hetlik(:,i,j,k,m)=bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j)*corr(:,i,j,k,m),W(j)*(1-corr(:,i,j,k,m)))+(corr(:,i,j,k,m)./0.5).*bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j)*(1-corr(:,i,j,k,m)),W(j)*corr(:,i,j,k,m))+inputParam.minLik;
                %hetlik(:,i,j,k)=bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j)*corr(:,i,j,k),W(j)*(1-corr(:,i,j,k)))+inputParam.minLik;
                hetlik(corr(:,i,j,k,m)==0,i,j,k,m)=inputParam.minLik;
                %hetlik(corr(:,i,j,k)==0.5,i,j,k)=bbinopdf_ln(D.MinorReadCount(corr(:,i,j,k)==0.5,j),D.TotalReadCount(corr(:,i,j,k)==0.5,j),W(j)*0.5,W(j)*0.5)+bbinopdf_ln(D.TotalReadCount(corr(:,i,j,k)==0.5,j)-D.MinorReadCount(corr(:,i,j,k)==0.5,j),D.TotalReadCount(corr(:,i,j,k)==0.5,j),W(j)*0.5,W(j)*0.5)+inputParam.minLik;
                %depthlik(:,i,j)=normpdf(log(E.TumorRD(:,j)+1),log(expReadCount(:,i,j)+1),0.6);
                pHetDetect(:,i,j,k,m)=binocdf(inputParam.minBCount,round(meanTumorRDexon(:,j)),corrSeg(:,i,j,k,m),'upper');
                hetCountLik(:,i,j,k,m)=binocdf(hetCount(:,i,k,m),dbCounts,pHetDetect(:,i,j,k,m)*pHet);
                hetExp(:,i,j,k,m)=dbCounts.*pHetDetect(:,i,j,k,m)*pHet;
                if ismember(j,inputParam.contamIdx)
                    hetlik(:,i,j,k,m)=NaN;
                    %hetCount(:,i,j,k)=NaN;
                end
                %segLik(:,i,j,k)=nansum([(hetCount./(hetExp(:,i,j,k)+1)).*getMeanInRegions([D.Chr D.Pos],log(hetlik(:,i,j,k))+log(priorCNAf(:,i,k))+log(priorMinAllele(:,i,k)),segsMerged) getMeanInRegions([E.Chr E.StartPos],log(depthlik(:,i,j,k))+log(priorCNA(:,i,k)),segsMerged) log(hetCountLik(:,i,j,k))],2);
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
            %         hetExpMax(cnaIdx==i & cnIdx==k,:,:)=hetExp(cnaIdx==i & cnIdx==k,i,:,k);
            %         for j=1:length(Tcell)
            %             hetlikMax(cnaIdx(idx)==i & cnIdx(idx)==k,j)=hetlik(cnaIdx(idx)==i & cnIdx(idx)==k,i,j,k);
            %             depthlikMax(cnaIdx(idxExon)==i & cnIdx(idxExon)==k,j)=depthlik(cnaIdx(idxExon)==i & cnIdx(idxExon)==k,i,j,k);
            %             hetCountLikMax(cnaIdx==i & cnIdx==k,j)=hetCountLik(cnaIdx==i & cnIdx==k,i,j,k);
            %         end
            %         priorCNAMax(cnaIdx(idxExon)==i & cnIdx(idxExon)==k,:)=priorCNA(cnaIdx(idxExon)==i & cnIdx(idxExon)==k,i,k);
            %         priorMinAlleleMax(cnaIdx(idx)==i & cnIdx(idx)==k,:)=priorMinAllele(cnaIdx(idx)==i & cnIdx(idx)==k,i,k);
            %         priorCNAfmax(cnaIdx(idx)==i & cnIdx(idx)==k,:)=priorCNAf(cnaIdx(idx)==i & cnIdx(idx)==k,i,k);
            %priorCNAfmax(cnaIdx(idx)==i & cnIdx(idx)==k,:)=betapdf(max(f(:,i)),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2)+inputParam.minLik;
        end
    end
end
% if inputParam.NormalSample>0
%     priorCNAfmax(cnaIdx(idx)==size(f,2),:)=germPriorMat(cnaIdx(idx)==size(f,2);
% end
%priorCNAfmax(NsegMax(idx)==2 & MsegMax(idx)==1)=NaN;

Fout=NaN(length(cnaIdx),length(Tcell));
for j=1:length(Tcell)
    Fout(:,j)=f(j,cnaIdx);
    %Wout(:,j)=W(j,cnaIdx);
end

nll=sum(segLikMax);
%nll=sum(-nansum(mean(log(hetlikMax)))-sum(mean(log(depthlikMax)))-sum((segsMerged(:,3)-segsMerged(:,2))'*log(hetCountLikMax)/sum(segsMerged(:,3)-segsMerged(:,2)))-mean(log(priorCNAMax))-mean(log(priorMinAlleleMax))-nanmean(log(priorCNAfmax)));
%%% missing somLik, priorF

return
