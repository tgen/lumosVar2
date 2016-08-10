function [N, M, Fout, log2FC] = callCNAmulti(hetPos,Tcell,exonRD,segsMerged,inputParam,param)
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
f=[zeros(length(Tcell),inputParam.numClones) ones(length(Tcell),1)];
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
f(tIdx,1:end-1)=reshape(param(2*length(Tcell)+1:end),[],inputParam.numClones);

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


%%% lookup copy number per position or exon
idx=getPosInRegions([D.Chr D.Pos], segsMerged);
Nmat=Nseg(idx,:);
Mmat=Mseg(idx,:);
idxExon=getPosInRegions([E.Chr E.StartPos],segsMerged);
NmatExon=Nseg(idxExon,:);

%%% find liklihoods of read counts and depth
for i=1:size(f,2)
    for j=1:length(Tcell)
        corr(:,i,j)=f(j,i).*Mmat(:,i)./Nmat(:,i)+(1-f(j,i))*0.5;
        corr(Nmat(:,i)==0,i,j)=0.5;
        hetlik(:,i,j)=bbinopdf_ln(D.MinorReadCount(:,j),D.TotalReadCount(:,j),W(j,:)*corr(:,i,j),W(j,:)*(1-corr(:,i,j)))+inputParam.minLik;
        expReadCount(:,i,j)=f(j,i)*E.NormalRD(:,j).*NmatExon(:,i)./CNAscale(j)+(1-f(j,i))*E.NormalRD(:,j)*2./CNAscale(j);
        depthlik(:,i,j)=poisspdf(round(E.TumorRD(:,j)),round(expReadCount(:,i,j)))+inputParam.minLik;
        segLik(:,i,j)=getMeanInRegions([D.Chr D.Pos],log(hetlik(:,i,j)),segsMerged)+getMeanInRegions([E.Chr E.StartPos],log(depthlik(:,i,j)),segsMerged);
    end
end


%%% find which clone contains most likely CNV per segment
[~,cnaIdx]=max(sum(segLik,3),[],2);
for i=1:size(f,2)
    %for j=1:length(Tcell)
        N(cnaIdx==i,:)=Nseg(cnaIdx==i,i);
        M(cnaIdx==i,:)=Mseg(cnaIdx==i,i);
    %end
end
for j=1:length(Tcell)
    Fout(:,j)=f(j,cnaIdx);
    %Wout(:,j)=W(j,cnaIdx);
end

return