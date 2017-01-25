function nll = nllCNAaddClone(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam,CNAscale,W,fOld,fNew)
%nllCNA - computes negative loglikliehood of copy number parameters
%
% Syntax: nll = nllCNA(dataHet,dataSom,exonRD,segs,inputParam,param)
%
% Inputs:
%   dataHet: data for germline heterozygous positions with columns:
%       1-'Chr',2-'Pos',3-'ControlRD',4-'TumorRD',5-'Bcount'
%   dataSom: data for somatic positions with columns:
%       1-'Chr',2-'Pos',3-'ControlRD',4-'TumorRD',5-'Bcount'
%   exonRD: matrix of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segs: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%   inputParam: structure with fields: cnaPrior, minAllelePrior, minLik,
%       alphaF, priorF
%   param: vector of paramaters of length 2*numClones+1
%       param(1) - copy number scaling constant
%       param(2:numClones+1)=W (controls width of allele frequency dist)
%       param(numClones+1:2*numClones+1)=f (sample fractions)
%   
% Outputs:
%   nll: negative log likelihood
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, fitCNA, callCNA

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

%%% read in inputs
D=table();
D.Chr=Tcell{1}.Chr(hetPos);
D.Pos=Tcell{1}.Pos(hetPos);
if(sum(somPos)>0)
    S=table();
    S.Chr=Tcell{1}.Chr(somPos);
    S.Pos=Tcell{1}.Pos(somPos);
end
E=table();
E.Chr=exonRD{1}(:,1);
E.StartPos=exonRD{1}(:,2);
E.EndPos=exonRD{1}(:,3);
for i=1:length(Tcell)
    D.ExpReadCount(:,i)=Tcell{i}.ControlRD(hetPos);
    D.TotalReadCount(:,i)=Tcell{i}.ReadDepthPass(hetPos);
    D.MinorReadCount(:,i)=Tcell{i}.BCountF(hetPos)+Tcell{i}.BCountR(hetPos);
    if(sum(somPos)>0)
        S.ExpReadCount(:,i)=Tcell{i}.ControlRD(somPos);
        S.TotalReadCount(:,i)=Tcell{i}.ReadDepthPass(somPos);
        S.MinorReadCount(:,i)=Tcell{i}.BCountF(somPos)+Tcell{i}.BCountR(somPos);
    end
    E.TumorRD(:,i)=exonRD{i}(:,4);
    E.NormalRD(:,i)=exonRD{i}(:,5);
end
if inputParam.NormalSample>0
    f=[zeros(length(Tcell),inputParam.numClones) ones(length(Tcell),1)];
    tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
    f(tIdx,1:end-1)=[fOld fNew(:)]./100;
else
    f=[fOld fNew(:)]./100;
end
    
%%% find means accross segments
for i=1:length(Tcell)
    meanTumorRDexon(:,i)=getMeanInRegions([E.Chr E.StartPos],E.TumorRD(:,i),segsMerged);
    meanNormalRDexon(:,i)=getMeanInRegions([E.Chr E.StartPos],E.NormalRD(:,i),segsMerged);
    meanTumorRD(:,i)=getMeanInRegions([D.Chr D.Pos],D.TotalReadCount(:,i),segsMerged);
    meanMinorRD(:,i)=getMeanInRegions([D.Chr D.Pos],D.MinorReadCount(:,i),segsMerged);
end

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


%%% lookup copy number for positions and exons
idx=getPosInRegions([D.Chr D.Pos], segsMerged);
Nmat=Nseg(idx,:,:);
Mmat=Mseg(idx,:,:);
idxExon=getPosInRegions([E.Chr E.StartPos],segsMerged);
NmatExon=Nseg(idxExon,:,:);

%%% find prior of copy number
priorCNA=nan(size(NmatExon));
for i=1:length(inputParam.cnaPrior)-1
    priorCNA(NmatExon==i-1)=inputParam.cnaPrior(i);
end
priorCNA(NmatExon>=length(inputParam.cnaPrior)-1)=inputParam.cnaPrior(end);
priorMinAllele=nan(size(Mmat));
for i=1:length(inputParam.minAllelePrior)-1
    priorMinAllele(Mmat==i-1)=inputParam.minAllelePrior(i);
end
priorMinAllele(Mmat>=length(inputParam.minAllelePrior)-1)=inputParam.minAllelePrior(end);

%%% find likelihoods of read counts and depth
pHet=inputParam.pvFreq*sum(E.EndPos-E.StartPos)./sum(segsMerged(:,3)-segsMerged(:,2));
priorCNAf=NaN(size(Nmat));
for i=1:size(f,2)
    for k=1:2
        if inputParam.NormalSample>0 && i==size(f,2)
            priorCNAf(:,i,k)=inputParam.priorGermCNV;
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
        priorCNAfmax(cnaIdx(idx)==i & cnIdx(idx)==k,:)=betapdf(max(f(:,i)),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2)+inputParam.minLik;
    end
end
if inputParam.NormalSample>0
    priorCNAfmax(cnaIdx(idx)==size(f,2),:)=inputParam.priorGermCNV;
end
priorCNAfmax(NsegMax(idx)==2 & MsegMax(idx)==1)=NaN;

%%% find expected allele frequency for somatic variants
for i=1:size(f,2)
    for j=1:length(Tcell)
        expAF(cnaIdx==i,i,j)=f(j,i)*(NsegMax(cnaIdx==i,:)-MsegMax(cnaIdx==i,:))./(f(j,i)*NsegMax(cnaIdx==i,:)+(1-f(j,i))*2);
        %subIdx=f(j,cnaIdx)+f(j,i)>1 & f(j,i)<f(j,cnaIdx);
        if sum(cnaIdx~=i)>0
            expAF(cnaIdx~=i,i,j)=f(j,i)./(f(j,cnaIdx(cnaIdx~=i))'.*NsegMax(cnaIdx~=i,:)+(1-f(j,cnaIdx(cnaIdx~=i))')*2);
         %   expAF(NsegMax==0 & cnaIdx~=i,i,j)=f(j,i)./2;
          %  expAF(subIdx,i,j)=f(j,i)*MsegMax(subIdx,:)./(f(j,cnaIdx(subIdx))'.*NsegMax(subIdx,:)+(1-f(j,cnaIdx(subIdx))')*2);
        end
        %[ones(sum(expAF(:,i,j)>1),1)*[i j f(j,i)] NsegMax(expAF(:,i,j)>1,:) MsegMax(expAF(:,i,j)>1,:) expAF(expAF(:,i,j)>1,i,j) cnaIdx(expAF(:,i,j)>1)]
        %[ones(sum(expAF(:,i,j)<0),1)*[i j f(j,i)] NsegMax(expAF(:,i,j)<0,:) MsegMax(expAF(:,i,j)<0,:) expAF(expAF(:,i,j)<0,i,j) cnaIdx(expAF(:,i,j)<0)]
    end
end
if(sum(somPos)>0)
    expAF(expAF>1)=1;
    expAF(expAF<0)=0;
    %%% find likelihood of somatic variant
    idxSom=getPosInRegions([S.Chr S.Pos], segsMerged);
    for i=1:size(f,2)
        for j=1:length(Tcell)
            alpha(:,i,j)=max(expAF(idxSom,i,j),inputParam.minLik)*W(j);
            beta(:,i,j)=(1-max(expAF(idxSom,i,j),inputParam.minLik))*W(j);
            cloneLik(:,i,j)=bbinopdf_ln(S.MinorReadCount(:,j),S.TotalReadCount(:,j),alpha(:,i,j),beta(:,i,j))+inputParam.minLik;
        end
    end
    cloneLik(isnan(cloneLik))=1;
    if (inputParam.NormalSample>0)
        cloneLik(:,end,:)=0;
        cloneLik(:,:,inputParam.NormalSample)=1;
    end
    [somLik,somIdx]=max(prod(cloneLik,3),[],2);
    %for j=1:length(Tcell)
    fMax=max(f,[],1);
    priorF=betapdf(fMax(somIdx),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2)+inputParam.minLik;
else
    somLik=1;
    priorF=1;
end
%end
%priorF=1;

%%% sum negative log likliehoods
%nll=sum((-sum(log(somLik))-sum(log(hetlikMax))-sum(log(depthlikMax))-sum(log(priorCNAMax))-sum(log(priorMinAlleleMax))-sum(log(priorF))-nansum(log(priorCNAf)))./(length(somLik)+length(hetlikMax)+length(depthlikMax)+length(priorCNAMax)+length(priorMinAlleleMax)+length(priorF)+sum(~isnan(priorCNAf))));
nll=sum(-sum(log(somLik)./(inputParam.priorSomaticSNV*sum(E.EndPos-E.StartPos)))-sum(mean(log(hetlikMax)))-sum(mean(log(depthlikMax)))-sum((segsMerged(:,3)-segsMerged(:,2))'*log(hetCountLikMax)/sum(segsMerged(:,3)-segsMerged(:,2)))-mean(log(priorCNAMax))-mean(log(priorMinAlleleMax))-mean(log(priorF))-nanmean(log(priorCNAfmax)));

return;
