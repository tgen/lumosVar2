function [F,pQual]=qualDiscrimAutoCut(T,E,homPos,inputParam)
%qualDiscrim - finds probability that positions are trusted or artifacts
%using quadratic discriminant analysis on a number of quality metrics
%
% Syntax:  [F,pTrust,pArtifact]=qualDiscrim(T,inputParam)
%
% Inputs:
%    T - table of position data with the following columns: 'Chr','Pos',
%       'ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ',
%       'AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ',
%       'BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount',
%       'ControlRD','PosMapQC','perReadPass','abFrac'
%    E - table of exon datawith the following columns: 'Chr','StartPos','EndPos',
%       'TumorRD','NormalRD', 'MapQC', 'perReadPass', 'abFrac'
%   inputParam - data structure with the following fields: ReadLength, 
%       minPerReadPASS, minABFrac, minPercentStrand, minMeanBQ, minMeanMQ, 
%       maxPMM, minSeqEndDist, maxStrandDiff, maxBQdiff, maxMQdiff, maxPMMdiff,
%       maxReadPosDiff, minPosQual, minExonQual, perPassReadReject,ABfracReject, 
%       perStrandReject, meanBQReject, meanMQReject, PMMReject,seqEndDistReject, 
%       strandDiffReject, BQdiffReject, MQdiffReject, PMMdiffReject, 
%       ReadPosDiffReject, inputParam.posQualReject, exonQualReject
%   
% Outputs:
%    F - table of quality statistics on positions same order as T
%    pTrust - posterior probability the call should be trusted
%    pArtifact - posterior probabilty the position has an artifact
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, preproccessTumorOnly

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------



%%%set quality metrics of positions with no B allele to NaN
T.BmeanBQ(T.BCountF+T.BCountR==0)=NaN;
T.BmeanMQ(T.BCountF+T.BCountR==0)=NaN;
T.BmeanPMM(T.BCountF+T.BCountR==0)=NaN;
T.BmeanReadPos(T.BCountF+T.BCountR==0)=NaN;

%%%determine whether homozygous positions are A or B
aIdx=homPos & T.ACountF+T.ACountR >= T.BCountF+T.BCountR;
bIdx=homPos & T.ACountF+T.ACountR < T.BCountF+T.BCountR;

%%%Create Quality Metrics Table
F=table();
F.TumorPerPassReads=-10*log10(1-T.ReadDepthPass./T.ReadDepth+1E-5);
F.TumorPerPassReads(T.ReadDepth==0)=0;
F.normalPerReadPass=-10*log10(1-T.perReadPass+1E-5);
F.normalPerReadPass(~isfinite(T.perReadPass))=-10*log10(0.5);
F.ABfrac=-10*log10(1-((T.ACountF+T.ACountR+T.BCountF+T.BCountR)./T.ReadDepthPass) + 1E-5);
F.normalABfrac=-10*log10(1-min(T.abFrac,1)+1E-5);
[F.minPerStrand,idx]=min([T.ACountF./(T.ACountF+T.ACountR) T.ACountR./(T.ACountF+T.ACountR) T.BCountF./(T.BCountF+T.BCountR) T.BCountR./(T.BCountF+T.BCountR)]+1E-5,[],2);
F.minPerStrand(aIdx)=min([T.ACountF(aIdx)./(T.ACountF(aIdx)+T.ACountR(aIdx)) T.ACountR(aIdx)./(T.ACountF(aIdx)+T.ACountR(aIdx)) 0.5*ones(sum(aIdx),1)],[],2);
F.minPerStrand(bIdx)=min([T.BCountF(bIdx)./(T.BCountF(bIdx)+T.BCountR(bIdx)) T.BCountR(bIdx)./(T.BCountF(bIdx)+T.BCountR(bIdx)) 0.5*ones(sum(bIdx),1)],[],2);
% pStrand=NaN(height(T),1);
% pStrand(idx==1)=arrayfun(@(x) hygecdf(T.ACountF(x),T.ReadDepthPass(x),round(T.ACountF(x)+T.BCountF(x)),T.ACountF(x)+T.ACountR(x)),find(idx==1));
% pStrand(idx==2)=arrayfun(@(x) hygecdf(T.ACountR(x),T.ReadDepthPass(x),round(T.ACountR(x)+T.BCountR(x)),T.ACountF(x)+T.ACountR(x)),find(idx==2));
% pStrand(idx==3)=arrayfun(@(x) hygecdf(T.BCountF(x),T.ReadDepthPass(x),round(T.ACountF(x)+T.BCountF(x)),T.BCountF(x)+T.BCountR(x)),find(idx==3));
% pStrand(idx==4)=arrayfun(@(x) hygecdf(T.BCountR(x),T.ReadDepthPass(x),round(T.ACountR(x)+T.BCountR(x)),T.BCountF(x)+T.BCountR(x)),find(idx==4));
F.minBQ=min([T.AmeanBQ T.BmeanBQ],[],2);
F.minBQ(aIdx)=T.AmeanBQ(aIdx);
F.minBQ(bIdx)=T.BmeanBQ(bIdx);
F.minMQ=min([T.AmeanMQ T.BmeanMQ],[],2);
F.minMQ(aIdx)=T.AmeanMQ(aIdx);
F.minMQ(bIdx)=T.BmeanMQ(bIdx);
F.maxPMM=-10*log10(max([T.AmeanPMM T.BmeanPMM]+1E-5,[],2));
F.maxPMM(aIdx)=-10*log10(T.AmeanPMM(aIdx)+1E-5);
F.maxPMM(bIdx)=-10*log10(T.BmeanPMM(bIdx)+1E-5);
inputParam.ReadLength=max([T.AmeanReadPos; T.BmeanReadPos]);
F.seqEndDist=min([T.AmeanReadPos inputParam.ReadLength-T.AmeanReadPos T.BmeanReadPos inputParam.ReadLength-T.BmeanReadPos],[],2);   
F.seqEndDist(aIdx)=min([T.AmeanReadPos(aIdx) inputParam.ReadLength-T.AmeanReadPos(aIdx)],[],2);   
F.seqEndDist(bIdx)=min([T.BmeanReadPos(bIdx) inputParam.ReadLength-T.BmeanReadPos(bIdx)],[],2);
F.strandDiff=-max(abs(T.ACountF./(T.ACountF+T.ACountR)-T.BCountF./(T.BCountF+T.BCountR)),0);
%F.strandDiff=10*log10(pStrand+1E-5);
F.strandDiff(homPos)=0;
F.BQdiff=-max(abs(T.AmeanBQ-T.BmeanBQ),0);
F.BQdiff(homPos)=0;
F.MQdiff=-max(abs(T.AmeanMQ-T.BmeanMQ),0);
F.MQdiff(homPos)=0;
F.PMMdiff=-10*log10(max(abs(T.AmeanPMM-T.BmeanPMM),1E-5));
F.PMMdiff(homPos)=0;
F.ReadPosDiff=-max(abs((T.AmeanReadPos-T.BmeanReadPos)),0);
F.ReadPosDiff(homPos)=0;
F.posMapQC=-10*log10(T.PosMapQC+1E-6);
F.posMapQC(T.PosMapQC<0)=60;
F.posMapQC(isnan(T.PosMapQC))=0;
idx=getPosInRegionSplit([T.Chr T.Pos],[E.Chr E.StartPos E.EndPos+1],inputParam.blockSize);
F.exonMapQC=zeros(height(F),1);
F.exonMapQC(~isnan(idx))=-10*log10(E.MapQC(idx(~isnan(idx)))+1E-6);
F.exonMapQC(E.MapQC(idx(~isnan(idx)))<0)=60;
F.exonMapQC(isnan(E.MapQC(idx(~isnan(idx)))))=0;
[F.ROmin,idx]=min([T.AmeanRO 1-T.AmeanRO T.BmeanRO 1-T.BmeanRO]+1E-5, [], 2);
F.ROmin(aIdx)=min([T.AmeanRO(aIdx) 1-T.AmeanRO(aIdx)],[],2);
F.ROmin(bIdx)=min([T.BmeanRO(bIdx) 1-T.BmeanRO(bIdx)],[],2);
T.ACountROF=T.AmeanRO.*(T.ACountF+T.ACountR);
T.ACountROR=(1-T.AmeanRO).*(T.ACountF+T.ACountR);
T.BCountROF=T.BmeanRO.*(T.BCountF+T.BCountR);
T.BCountROR=(1-T.BmeanRO).*(T.BCountF+T.BCountR);
% pRO(idx==1)=arrayfun(@(x) hygecdf(T.ACountROF(x),T.ReadDepthPass(x),round(T.ACountROF(x)+T.BCountROF(x)),round(T.ACountROF(x)+T.ACountROR(x))),find(idx==1));
% pRO(idx==2)=arrayfun(@(x) hygecdf(T.ACountROR(x),T.ReadDepthPass(x),round(T.ACountROR(x)+T.BCountROR(x)),round(T.ACountROF(x)+T.ACountROR(x))),find(idx==2));
% pRO(idx==3)=arrayfun(@(x) hygecdf(T.BCountROF(x),T.ReadDepthPass(x),round(T.ACountROF(x)+T.BCountROF(x)),round(T.BCountROF(x)+T.BCountROR(x))),find(idx==3));
% pRO(idx==4)=arrayfun(@(x) hygecdf(T.BCountROR(x),T.ReadDepthPass(x),round(T.ACountROR(x)+T.BCountROR(x)),round(T.BCountROF(x)+T.BCountROR(x))),find(idx==4));
% F.ROdiff=10*log10(pRO'+1E-5);
F.ROdiff=-max(abs(T.AmeanRO-T.BmeanRO),0);
F.ROdiff(homPos)=0;
F.ISmax=-max([T.AmeanIS T.BmeanIS], [], 2);
F.ISmax(aIdx)=-T.AmeanIS(aIdx);
F.ISmax(bIdx)=-T.BmeanIS(bIdx);
F.ISmin=min([T.AmeanIS T.BmeanIS], [], 2);
F.ISmin(aIdx)=T.AmeanIS(aIdx);
F.ISmin(bIdx)=T.BmeanIS(bIdx);
F.ISdiff=-max(abs(T.AmeanIS - T.BmeanIS), 0);
F.ISdiff(homPos)=0;
F.SCmax=-max([T.AmeanSC T.BmeanSC], [], 2);
F.SCmax(aIdx)=-T.AmeanSC(aIdx);
F.SCmax(bIdx)=-T.BmeanSC(bIdx);
F.SCdiff=-max(abs(T.AmeanSC - T.BmeanSC), 0);
F.SCdiff(homPos)=0;

F.Properties.VariableDescriptions={'Percent of Reads in Tumor Passing Quality Thresh', ...
    'Percent of Reads in Normals Passing Quality Thresh', ...
    'Fraction of QC Reads in Tumor Supporting A or B Allele', ...
    'Fraction of QC Reads in Normals Supporting A or B Allele', ...
    'Minimum Percentage of Reads from forward or reverse strand supporting A or B allele', ...
    'Minimum of average BQ for reads supporing A or B allele', ...
    'Minimum of average MQ for reads supporting A or B allele', ...
    'Maximum percentage of mismatches in reads supporting A or B allele', ...
    'Minimum average distance to end of sequence of reads supporting A or B allele', ...
    'Phred scaled Fisher test for strand bias', ...
    'Difference in average BQ between bases supporting A and B allele', ...
    'Difference in average MQ between reads supporting A and B allele', ...
    'Difference in average percentage of mismatches in reads supporting A and B allele', ...
    'Difference in average position in sequence between A and B allele', ...
    'Phred scaled postion quality score from Normals', ...
    'mean position quality score in exon','','','','','','',''};

[~,lb_good_nh]=isoutlier(F{~homPos &F.posMapQC>30,:},'mean','ThresholdFactor',inputParam.outThreshGood);
%lb_good_nh=nanmedian(F{~homPos &F.posMapQC>30,:})-inputParam.outThreshGood*mad(F{~homPos &F.posMapQC>30,:},1)*1.4862;
[~,lb_good_h]=isoutlier(F{homPos &F.posMapQC>30,:},'mean','ThresholdFactor',inputParam.outThreshGood);
%lb_good_h=nanmedian(F{homPos &F.posMapQC>30,:})-inputParam.outThreshGood*mad(F{homPos &F.posMapQC>30,:},1)*1.4862;


goodPos=zeros(size(F));
goodPos(~homPos,:)=F{~homPos,:}>=ones(sum(~homPos),1)*lb_good_nh;
goodPos(homPos,:)=F{homPos,:}>=ones(sum(homPos),1)*lb_good_h;
%sum(goodPos(~homPos,:))./sum(~homPos)
%histogram(sum(goodPos,2))

[~,lb_bad_nh]=isoutlier(F{~homPos &F.posMapQC>30,:},'mean','ThresholdFactor',inputParam.outThreshBad);
[~,lb_bad_h]=isoutlier(F{homPos &F.posMapQC>30,:},'mean','ThresholdFactor',inputParam.outThreshBad);
%lb_bad_nh=nanmedian(F{~homPos &F.posMapQC>30,:})-inputParam.outThreshBad*mad(F{~homPos &F.posMapQC>30,:},1)*1.4862;
%lb_bad_h=nanmedian(F{homPos &F.posMapQC>30,:})-inputParam.outThreshBad*mad(F{homPos &F.posMapQC>30,:},1)*1.4862;


badPos=zeros(size(F));
badPos(~homPos,:)=F{~homPos,:}<ones(sum(~homPos),1)*lb_bad_nh;
badPos(homPos,:)=F{homPos,:}<ones(sum(homPos),1)*lb_bad_h;

mdl=fitcdiscr([F{~homPos,:}; F{~homPos,:}],[ones(sum(~homPos),1); zeros(sum(~homPos),1)],'DiscrimType','pseudoquadratic','Weights',[1./(size(goodPos,2)+1-sum(goodPos(~homPos,:),2)); log(sum(badPos(~homPos,:),2)+1)]);
[~,pQual]=predict(mdl,F{:,:});
pQual(T.ReadDepthPass==0,:)=0.5;

F.goodPosSum=sum(goodPos,2);
F.badPosSum=sum(badPos,2);

F.Properties.VariableDescriptions(24)={'Number of PASS criteria met'};
F.Properties.VariableDescriptions(25)={'Number of REJECT criteria met'};
return;
