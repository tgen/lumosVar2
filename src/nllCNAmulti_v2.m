function nll = nllCNAmulti_v2(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam,param)
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
S=table();
S.Chr=Tcell{1}.Chr(somPos);
S.Pos=Tcell{1}.Pos(somPos);
E=table();
E.Chr=exonRD{1}(:,1);
E.StartPos=exonRD{1}(:,2);
E.EndPos=exonRD{1}(:,3);
for i=1:length(Tcell)
    D.ExpReadCount(:,i)=Tcell{i}.ControlRD(hetPos);
    D.TotalReadCount(:,i)=Tcell{i}.ReadDepthPass(hetPos);
    D.MinorReadCount(:,i)=Tcell{i}.BCountF(hetPos)+Tcell{i}.BCountR(hetPos);
    S.ExpReadCount(:,i)=Tcell{i}.ControlRD(somPos);
    S.TotalReadCount(:,i)=Tcell{i}.ReadDepthPass(somPos);
    S.MinorReadCount(:,i)=Tcell{i}.BCountF(somPos)+Tcell{i}.BCountR(somPos);
    E.TumorRD(:,i)=exonRD{i}(:,4);
    E.NormalRD(:,i)=exonRD{i}(:,5);
end
CNAscale=param(1:length(Tcell))./100;
W=param(length(Tcell)+1:2*length(Tcell));
f=[zeros(length(Tcell),inputParam.numClones) ones(length(Tcell),1)];
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
f(tIdx,1:end-1)=reshape(param(2*length(Tcell)+1:end)./100,[],inputParam.numClones);

%%% find means accross segments
for i=1:length(Tcell)
    meanTumorRDexon(:,i)=getMeanInRegions([E.Chr E.StartPos],E.TumorRD(:,i),segsMerged);
    meanNormalRDexon(:,i)=getMeanInRegions([E.Chr E.StartPos],E.NormalRD(:,i),segsMerged);
    meanTumorRD(:,i)=getMeanInRegions([D.Chr D.Pos],D.TotalReadCount(:,i),segsMerged);
    meanMinorRD(:,i)=getMeanInRegions([D.Chr D.Pos],D.MinorReadCount(:,i),segsMerged);
end

%%% find most likely copy number and minor for each segment and clone
for i=1:size(f,2)
    for j=1:length(Tcell)
        if(f(j,i)==0)
            NsegSample(:,i,j)=2*ones(size(segsMerged,1),1);
            MsegSample(:,i,j)=ones(size(segsMerged,1),1);
        else
            NsegSample(:,i,j)=max((CNAscale(j)*(meanTumorRDexon(:,j)./(meanNormalRDexon(:,j)))-2*(1-f(j,i)))/f(j,i),0);
            MsegSample(isfinite(meanMinorRD(:,j)),i,j)=(NsegSample(isfinite(meanMinorRD(:,j)),i,j)/f(j,i)).*((meanMinorRD(isfinite(meanMinorRD(:,j)),j)./(meanTumorRD(isfinite(meanMinorRD(:,j)),j)))-0.5*(1-f(j,i)));
            MsegSample(~isfinite(meanMinorRD(:,j)),i,j)=min(NsegSample(~isfinite(meanMinorRD(:,j)),i,j)-1,1);
        end
    end
end
NsegSample(NsegSample<0)=0;
MsegSample(MsegSample<0)=0;
for i=1:size(f,2)
    Nseg(:,i)=round(squeeze(NsegSample(:,i,:))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
    Mseg(:,i)=round(squeeze(MsegSample(:,i,:))*f(:,i)./(ones(size(segsMerged,1),size(f,1))*f(:,i)));
end


%%% lookup copy number for positions and exons
idx=getPosInRegions([D.Chr D.Pos], segsMerged);
Nmat=Nseg(idx,:);
Mmat=Mseg(idx,:);
idxExon=getPosInRegions([E.Chr E.StartPos],segsMerged);
NmatExon=Nseg(idxExon,:);

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
for i=1:size(f,2)
    for j=1:length(Tcell)
        corr(:,i,j)=f(j,i).*Mmat(:,i)./Nmat(:,i)+(1-f(j,i))*0.5;
        corr(Nmat(:,i)==0,i,j)=0.5;
        corr(corr(:,i,j)<0,i,j)=0;
        corr(corr(:,i,j)>1,i,j)=1;
        if min(corr(:,i,j))<0 || max(corr(:,i,j))>1
            idx=corr(:,i,j)<0 | corr(:,i,j)>1;
            size(idx)
            [min(corr(:,i,j)) max(corr(:,i,j))]
            f(j,i)
%             %corr(idx,i,j) 
%             %Nmat(idx,i) 
%             %Mmat(idx,i)
         end
        hetlik(:,i,j)=bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j)*corr(:,i,j),W(j)*(1-corr(:,i,j)))+bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j)*(1-corr(:,i,j)),W(j)*corr(:,i,j))+inputParam.minLik;
        hetlik(corr(:,i,j)==0,i,j)=inputParam.minLik;
        expReadCount(:,i,j)=f(j,i)*E.NormalRD(:,j).*NmatExon(:,i)./CNAscale(j)+(1-f(j,i))*E.NormalRD(:,j)*2./CNAscale(j);
        depthlik(:,i,j)=poisspdf(round(E.TumorRD(:,j)),round(expReadCount(:,i,j)))+inputParam.minLik;
        segLik(:,i,j)=getMeanInRegions([D.Chr D.Pos],log(hetlik(:,i,j)),segsMerged)+getMeanInRegions([E.Chr E.StartPos],log(depthlik(:,i,j)),segsMerged);
    end
end

%%% find which clone contains most likely CNV per segment
[m,cnaIdx]=max(sum(segLik,3),[],2);
for i=1:size(f,2)
    NsegMax(cnaIdx==i,:)=Nseg(cnaIdx==i,i);
    MsegMax(cnaIdx==i,:)=Mseg(cnaIdx==i,i);
    for j=1:length(Tcell)   
        hetlikMax(cnaIdx(idx)==i,j)=hetlik(cnaIdx(idx)==i,i,j);
        depthlikMax(cnaIdx(idxExon)==i,j)=depthlik(cnaIdx(idxExon)==i,i,j); 
    end
    priorCNAMax(cnaIdx(idxExon)==i,:)=priorCNA(cnaIdx(idxExon)==i,i);
    priorMinAlleleMax(cnaIdx(idx)==i,:)=priorMinAllele(cnaIdx(idx)==i,i);
    priorCNAf(cnaIdx(idx)==i,:)=betapdf(max(f(:,i)),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2)+inputParam.minLik;
end
priorCNAf(NsegMax==2 & MsegMax==1)=NaN;
%priorCNAf=1;

%%% find expected allele frequency for somatic variants
for i=1:size(f,2)
    for j=1:length(Tcell)
        expAF(cnaIdx==i,i,j)=f(j,i)*(NsegMax(cnaIdx==i,:)-MsegMax(cnaIdx==i,:))./(f(j,i)*NsegMax(cnaIdx==i,:)+(1-f(j,i))*2);
        subIdx=f(j,cnaIdx)+f(j,i)>1 & f(j,i)<f(j,cnaIdx);
        if sum(cnaIdx~=i)>0
            expAF(cnaIdx~=i,i,j)=f(j,i)./(f(j,cnaIdx(cnaIdx~=i))'.*NsegMax(cnaIdx~=i,:)+(1-f(j,cnaIdx(cnaIdx~=i))')*2);
            expAF(subIdx,i,j)=f(j,i)*MsegMax(subIdx,j)./(f(j,cnaIdx(subIdx)).*NsegMax(subIdx,j)+(1-f(j,cnaIdx(subIdx)))*2);
            %expAF(NsegMax==0 & cnaIdx~=i,i,j)=min([(1-f(j,cnaIdx(cnaIdx~=i))'); f(j,i)])./2;
        end
        %[ones(sum(expAF(:,i,j)>1),1)*[i j f(j,i)] NsegMax(expAF(:,i,j)>1,j) MsegMax(expAF(:,i,j)>1,j) expAF(expAF(:,i,j)>1,i,j) cnaIdx(expAF(:,i,j)>1)]
        %[ones(sum(expAF(:,i,j)<0),1)*[i j f(j,i)] NsegMax(expAF(:,i,j)<0,j) MsegMax(expAF(:,i,j)<0,j) expAF(expAF(:,i,j)<0,i,j) cnaIdx(expAF(:,i,j)<0)]
    end
end

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
fMax=max(f);
priorF=betapdf(fMax(somIdx),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2)+inputParam.minLik;
%end
%priorF=1;

%%% sum negative log likliehoods
%nll=sum((-sum(log(somLik))-sum(log(hetlikMax))-sum(log(depthlikMax))-sum(log(priorCNAMax))-sum(log(priorMinAlleleMax))-sum(log(priorF))-nansum(log(priorCNAf)))./(length(somLik)+length(hetlikMax)+length(depthlikMax)+length(priorCNAMax)+length(priorMinAlleleMax)+length(priorF)+sum(~isnan(priorCNAf))));

nll=sum(-sum(log(somLik)./(inputParam.priorSomatic*sum(E.EndPos-E.StartPos)))-mean(log(hetlikMax))-mean(log(depthlikMax))-mean(log(priorCNAMax))-mean(log(priorMinAlleleMax))-mean(log(priorF))-nanmean(log(priorCNAf)));

return;