function nll = nllCNA(dataHet,dataSom,exonRD,segs,inputParam,param)
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
colHeaders={'Chr' 'Pos' 'ExpReadCount' 'TotalReadCount' 'MinorReadCount'};
D=array2table(dataHet,'VariableNames',colHeaders);
if ~isempty(dataSom)
    S=array2table(dataSom,'VariableNames',colHeaders);
end
E=array2table(exonRD,'VariableNames',{'Chr','StartPos','EndPos','TumorRD','NormalRD', 'MapQC', 'perReadPass', 'abFrac'});
CNAscale=param(1);
W(1,:)=param(2:(length(param)-1)./2+1);
f(1,:)=param((length(param)-1)./2+2:end);

%%% find means accross segments
meanTumorRDexon=getMeanInRegions([E.Chr E.StartPos],E.TumorRD,segs);
meanNormalRDexon=getMeanInRegions([E.Chr E.StartPos],E.NormalRD,segs);
meanTumorRD=getMeanInRegions([D.Chr D.Pos],D.TotalReadCount,segs);
meanMinorRD=getMeanInRegions([D.Chr D.Pos],D.MinorReadCount,segs);

%%% find most likely copy number and minor for each segment and clone
for i=1:length(f)
    Nseg(:,i)=max(round((CNAscale*(meanTumorRDexon./(meanNormalRDexon))-2*(1-f(i)))/f(i)),0);
    Mseg(isfinite(meanMinorRD),i)=round((Nseg(isfinite(meanMinorRD),i)/f(i)).*((meanMinorRD(isfinite(meanMinorRD))./(meanTumorRD(isfinite(meanMinorRD))))-0.5*(1-f(i))));
    Mseg(~isfinite(meanMinorRD),i)=min(Nseg(~isfinite(meanMinorRD),i)-1,1);
end
Nseg(Nseg<0)=0;
Mseg(Mseg<0)=0;

%%% lookup copy number for positions and exons
idx=getPosInRegions([D.Chr D.Pos], segs);
Nmat=Nseg(idx,:);
Mmat=Mseg(idx,:);
idxExon=getPosInRegions([E.Chr E.StartPos],segs);
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
for i=1:length(f)
    corr(:,i)=f(i).*Mmat(:,i)./Nmat(:,i)+(1-f(i))*0.5;
    corr(Nmat(:,i)==0,i)=0.5;
    hetlik(:,i)=bbinopdf_ln(D.MinorReadCount,D.TotalReadCount,W(i)*corr(:,i),W(i)*(1-corr(:,i)))+bbinopdf_ln(D.MinorReadCount,D.TotalReadCount,W(i)*(1-corr(:,i)),W(i)*corr(:,i))+inputParam.minLik;
    expReadCount(:,i)=f(i)*E.NormalRD.*NmatExon(:,i)./CNAscale+(1-f(i))*E.NormalRD*2./CNAscale;
    depthlik(:,i)=poisspdf(round(E.TumorRD),round(expReadCount(:,i)))+inputParam.minLik;
    segLik(:,i)=getMeanInRegions([D.Chr D.Pos],log(hetlik(:,i)),segs)+getMeanInRegions([E.Chr E.StartPos],log(depthlik(:,i)),segs);
end

%%% find which clone contains most likely CNV per segment
[m,cnaIdx]=max(segLik,[],2);
for i=1:length(f)
    NsegMax(cnaIdx==i,:)=Nseg(cnaIdx==i,i);
    MsegMax(cnaIdx==i,:)=Mseg(cnaIdx==i,i);
    hetlikMax(cnaIdx(idx)==i,:)=hetlik(cnaIdx(idx)==i,i);  
    depthlikMax(cnaIdx(idxExon)==i,:)=depthlik(cnaIdx(idxExon)==i,i);   
    priorCNAMax(cnaIdx(idxExon)==i,:)=priorCNA(cnaIdx(idxExon)==i,i);
    priorMinAlleleMax(cnaIdx(idx)==i,:)=priorMinAllele(cnaIdx(idx)==i,i);
    priorCNAf(cnaIdx(idx)==i,:)=betapdf(f(i),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2);
end
priorCNAf(NsegMax==2 & MsegMax==1)=NaN;

%%% find expected allele frequency for somatic variants
for i=1:length(f)
    expAF(cnaIdx==i,i)=f(i)*(NsegMax(cnaIdx==i,:)-MsegMax(cnaIdx==i,:))./(f(i)*NsegMax(cnaIdx==i,:)+(1-f(i))*2);
    subIdx=f(cnaIdx)+f(i)>1 & f(i)<f(cnaIdx);
    if sum(cnaIdx~=i)>0
        expAF(cnaIdx~=i,i)=f(i)./(f(cnaIdx(cnaIdx~=i))'.*NsegMax(cnaIdx~=i,:)+(1-f(cnaIdx(cnaIdx~=i))')*2);
        expAF(subIdx,i)=f(i)*MsegMax(subIdx)./(f(cnaIdx(subIdx)).*NsegMax(subIdx)+(1-f(cnaIdx(subIdx)))*2);
        %expAF(NsegMax==0 & cnaIdx~=i,i,j)=min([(1-f(j,cnaIdx(cnaIdx~=i))'); f(j,i)])./2;
    end
end

%%% find likelihood of somatic variant
if isempty(dataSom)
    somLik=1;
    priorF=1;
else
    idxSom=getPosInRegions([S.Chr S.Pos], segs);
    for i=1:length(f)
        alpha(:,i)=max(expAF(idxSom,i),inputParam.minLik)*W(i);
        beta(:,i)=(1-max(expAF(idxSom,i),inputParam.minLik))*W(i);
        cloneLik(:,i)=bbinopdf_ln(S.MinorReadCount,S.TotalReadCount,alpha(:,i),beta(:,i))+inputParam.minLik;
    end
    %cloneLik(isnan(cloneLik))=1;
    [somLik,somIdx]=max(cloneLik,[],2);
    priorF=betapdf(f(somIdx),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2);
end

%%% sum negative log likliehoods
nll=(-sum(log(somLik))-sum(log(hetlikMax))-sum(log(depthlikMax))-sum(log(priorCNAMax))-sum(log(priorMinAlleleMax))-sum(log(priorF))-nansum(log(priorCNAf)))./(length(somLik)+length(hetlikMax)+length(depthlikMax)+length(priorCNAMax)+length(priorMinAlleleMax)+length(priorF)+sum(~isnan(priorCNAf)));

return;