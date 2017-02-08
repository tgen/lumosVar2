function [NsegMax, MsegMax, Fout, log2FC, cnaIdx, nll] = callCNAmulti(hetPos,Tcell,exonRD,segsMerged,inputParam,param)
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
    D.TotalReadCount(:,i)=Tcell{i}.ReadDepthPass(hetPos);
    D.MinorReadCount(:,i)=Tcell{i}.BCountF(hetPos)+Tcell{i}.BCountR(hetPos);
    E.TumorRD(:,i)=exonRD{i}(:,4);
    E.NormalRD(:,i)=exonRD{i}(:,5);
end
CNAscale=param(1:length(Tcell));
W=param(length(Tcell)+1:2*length(Tcell));
if inputParam.NormalSample>0
    f=[zeros(length(Tcell),inputParam.numClones) ones(length(Tcell),1)];
    tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
    f(tIdx,1:end-1)=reshape(param(2*length(Tcell)+1:end),[],inputParam.numClones);
else
    f=[reshape(param(2*length(Tcell)+1:end),[],inputParam.numClones) ones(length(Tcell),1)];
end

%%% find means accross segments
for i=1:length(Tcell)
    meanTumorRDexon(:,i)=getMeanInRegions([E.Chr E.StartPos],E.TumorRD(:,i),segsMerged);
    meanNormalRDexon(:,i)=getMeanInRegions([E.Chr E.StartPos],E.NormalRD(:,i),segsMerged);
    meanTumorRD(:,i)=getMeanInRegions([D.Chr D.Pos],D.TotalReadCount(:,i),segsMerged);
    meanMinorRD(:,i)=getMeanInRegions([D.Chr D.Pos],D.MinorReadCount(:,i),segsMerged);
end

%%% calculate log2FC
for i=1:length(Tcell)
    log2FC(:,i)=log2((CNAscale(i)./2).*meanTumorRDexon(:,i)./(meanNormalRDexon(:,i)));
end

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
for i=1:size(f,2)
    Nseg(:,i,1)=floor(squeeze(NsegSample(:,i,:))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
    Nseg(:,i,2)=Nseg(:,i,1)+1;
    %Mseg(:,i)=round(squeeze(MsegSample(:,i,:))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
end
for i=1:size(f,2)
    for j=1:length(Tcell)
        if(f(j,i)==0)
            %NsegSample(:,i,j)=2*ones(size(segsMerged,1),1);
            MsegSample(:,i,j,1)=ones(size(segsMerged,1),1);
            MsegSample(:,i,j,2)=ones(size(segsMerged,1),1);
        else
            %NsegSample(:,i,j)=max((CNAscale(j)*(meanTumorRDexon(:,j)./(meanNormalRDexon(:,j)))-2*(1-f(j,i)))/f(j,i),0);
            MsegSample(isfinite(meanMinorRD(:,j)),i,j,1)=(Nseg(isfinite(meanMinorRD(:,j)),i,1)/f(j,i)).*((meanMinorRD(isfinite(meanMinorRD(:,j)),j)./(meanTumorRD(isfinite(meanMinorRD(:,j)),j)))-0.5*(1-f(j,i)));
            MsegSample(isfinite(meanMinorRD(:,j)),i,j,2)=(Nseg(isfinite(meanMinorRD(:,j)),i,2)/f(j,i)).*((meanMinorRD(isfinite(meanMinorRD(:,j)),j)./(meanTumorRD(isfinite(meanMinorRD(:,j)),j)))-0.5*(1-f(j,i)));
            MsegSample(~isfinite(meanMinorRD(:,j)),i,j,1)=min(Nseg(~isfinite(meanMinorRD(:,j)),i,1)-1,1);
            MsegSample(~isfinite(meanMinorRD(:,j)),i,j,2)=min(Nseg(~isfinite(meanMinorRD(:,j)),i,2)-1,1);
        end
    end
end
MsegSample(MsegSample<0)=0;
for i=1:size(f,2)
    %Nseg(:,i,1)=floor(squeeze(NsegSample(:,i,:))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
    %Nseg(:,i,2)=Nseg(:,i,1)+1;
    Mseg(:,i,1)=round(squeeze(MsegSample(:,i,:,1))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
    Mseg(:,i,2)=round(squeeze(MsegSample(:,i,:,2))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
end


germPrior=inputParam.priorGermCNV*ones(size(Nseg,1),2);
[lia,locb]=ismember(segsMerged(:,1),inputParam.autosomes);
germPrior(lia & Nseg(:,end,1)==2 & Mseg(:,end,1)==1,1)=1-inputParam.priorGermCNV;
germPrior(lia & Nseg(:,end,2)==2 & Mseg(:,end,2)==1,2)=1-inputParam.priorGermCNV;
sexChr=regexp(inputParam.sexChr,',','split');
chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr'];
for i=1:length(sexChr)
    chrIdx=find(strcmp(sexChr(i),chrList));
    lia=ismember(Nseg(:,end,1),inputParam.(sexChr{i}));
    germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,1)-Mseg(:,end,1)<=1,1)=0.5;
    lia=ismember(Nseg(:,end,2),inputParam.(sexChr{i}));
    germPrior(chrIdx==segsMerged(:,1) & lia & Nseg(:,end,2)-Mseg(:,end,2)<=1,2)=0.5;
end
    

%%% lookup copy number for positions and exons
idx=getPosInRegions([D.Chr D.Pos], segsMerged);
Nmat=Nseg(idx,:,:);
Mmat=Mseg(idx,:,:);
germPriorMat=germPrior(idx,:);
idxExon=getPosInRegions([E.Chr E.StartPos],segsMerged);
NmatExon=Nseg(idxExon,:,:);
germPriorExon=germPrior(idxExon,:);


%%% find prior of copy number
priorCNA=nan(size(NmatExon));
for i=1:length(inputParam.cnaPrior)-1
    priorCNA(NmatExon==i-1)=inputParam.cnaPrior(i);
end
priorCNA(NmatExon>=length(inputParam.cnaPrior)-1)=inputParam.cnaPrior(end);
priorCNA(:,end,:)=germPriorExon;
priorMinAllele=nan(size(Mmat));
for i=1:length(inputParam.minAllelePrior)-1
    priorMinAllele(Mmat==i-1)=inputParam.minAllelePrior(i);
end
priorMinAllele(Mmat>=length(inputParam.minAllelePrior)-1)=inputParam.minAllelePrior(end);
priorMinAllele(:,end,:)=germPriorMat;

%%% find likelihoods of read counts and depth
pHet=inputParam.pvFreq*sum(E.EndPos-E.StartPos)./sum(segsMerged(:,3)-segsMerged(:,2));
priorCNAf=NaN(size(Nmat));
for i=1:size(f,2)
    for k=1:2
        if i==size(f,2)
            priorCNAf(:,i,k)=germPriorMat(:,k);
        else
            priorCNAf(:,i,k)=betapdf(max(f(:,i)),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2)+inputParam.minLik;
        end
        priorCNAf(Nmat(:,i,k)==2 & Mmat(:,i,k)==1,i,k)=betapdf(inputParam.priorF,inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2);
    end
    for j=1:length(Tcell)
        for k=1:2
            corr(:,i,j,k)=f(j,i).*Mmat(:,i,k)./Nmat(:,i,k)+(1-f(j,i))*0.5;
            corr(Nmat(:,i,k)==0,i,j,k)=0.5;
            corr(corr(:,i,j,k)<0,i,j,k)=0;
            corr(corr(:,i,j,k)>1,i,j,k)=1;
            corrSeg(:,i,j,k)=f(j,i).*Mseg(:,i,k)./Nseg(:,i,k)+(1-f(j,i))*0.5;
            corrSeg(Nseg(:,i,k)==0,i,j,k)=0.5;
            corrSeg(corrSeg(:,i,j,k)<0,i,j,k)=0;
            corrSeg(corrSeg(:,i,j,k)>1,i,j,k)=1;
            hetlik(:,i,j,k)=bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j)*corr(:,i,j,k),W(j)*(1-corr(:,i,j,k)))+bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j)*(1-corr(:,i,j,k)),W(j)*corr(:,i,j,k))+inputParam.minLik;
            hetlik(corr(:,i,j,k)==0,i,j,k)=inputParam.minLik;
            expReadCount(:,i,j,k)=f(j,i)*E.NormalRD(:,j).*NmatExon(:,i,k)./CNAscale(j)+(1-f(j,i))*E.NormalRD(:,j)*2./CNAscale(j);
            depthlik(:,i,j,k)=poisspdf(round(E.TumorRD(:,j)),round(expReadCount(:,i,j,k)))+inputParam.minLik;
            %depthlik(:,i,j)=normpdf(log(E.TumorRD(:,j)+1),log(expReadCount(:,i,j)+1),0.6);
            pHetDetect(:,i,j,k)=binocdf(inputParam.minBCount,round(meanTumorRDexon(:,j)),corrSeg(:,i,j,k),'upper');
            hetCountLik(:,i,j,k)=binocdf(hist(idx,1:size(segsMerged,1))',round(segsMerged(:,3)-segsMerged(:,2)),pHetDetect(:,i,j,k)*pHet);
            segLik(:,i,j,k)=nansum([getMeanInRegions([D.Chr D.Pos],log(hetlik(:,i,j,k))+log(priorCNAf(:,i,k))+log(priorMinAllele(:,i,k)),segsMerged) getMeanInRegions([E.Chr E.StartPos],log(depthlik(:,i,j,k))+log(priorCNA(:,i,k)),segsMerged) log(hetCountLik(:,i,j,k))],2);
        end
    end
end

%%% find which clone contains most likely CNV per segment
[segLikMax,maxIdx]=max(sum(segLik,3),[],4);
[m,cnaIdx]=max(segLikMax,[],2);
for i=1:size(f,2)
    cnIdx(cnaIdx==i,:)=maxIdx(cnaIdx==i,i);
end

for i=1:size(f,2)
    for k=1:2
        NsegMax(cnaIdx==i & cnIdx==k,:)=Nseg(cnaIdx==i & cnIdx==k,i,k);
        MsegMax(cnaIdx==i & cnIdx==k,:)=Mseg(cnaIdx==i & cnIdx==k,i,k);
        for j=1:length(Tcell)
            hetlikMax(cnaIdx(idx)==i & cnIdx(idx)==k,j)=hetlik(cnaIdx(idx)==i & cnIdx(idx)==k,i,j,k);
            depthlikMax(cnaIdx(idxExon)==i & cnIdx(idxExon)==k,j)=depthlik(cnaIdx(idxExon)==i & cnIdx(idxExon)==k,i,j,k);
            hetCountLikMax(cnaIdx==i & cnIdx==k,j)=hetCountLik(cnaIdx==i & cnIdx==k,i,j,k);
        end
        priorCNAMax(cnaIdx(idxExon)==i & cnIdx(idxExon)==k,:)=priorCNA(cnaIdx(idxExon)==i & cnIdx(idxExon)==k,i,k);
        priorMinAlleleMax(cnaIdx(idx)==i & cnIdx(idx)==k,:)=priorMinAllele(cnaIdx(idx)==i & cnIdx(idx)==k,i,k);
        priorCNAfmax(cnaIdx(idx)==i & cnIdx(idx)==k,:)=priorCNAf(cnaIdx(idx)==i & cnIdx(idx)==k,i,k);
        %priorCNAfmax(cnaIdx(idx)==i & cnIdx(idx)==k,:)=betapdf(max(f(:,i)),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2)+inputParam.minLik;
    end
end
% if inputParam.NormalSample>0
%     priorCNAfmax(cnaIdx(idx)==size(f,2),:)=germPriorMat(cnaIdx(idx)==size(f,2);
% end
priorCNAfmax(NsegMax(idx)==2 & MsegMax(idx)==1)=NaN;

for j=1:length(Tcell)
    Fout(:,j)=f(j,cnaIdx);
    %Wout(:,j)=W(j,cnaIdx);
end

nll=sum(-sum(mean(log(hetlikMax)))-sum(mean(log(depthlikMax)))-sum((segsMerged(:,3)-segsMerged(:,2))'*log(hetCountLikMax)/sum(segsMerged(:,3)-segsMerged(:,2)))-mean(log(priorCNAMax))-mean(log(priorMinAlleleMax))-nanmean(log(priorCNAfmax)));
%%% missing somLik, priorF

return