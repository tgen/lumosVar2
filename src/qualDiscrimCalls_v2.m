function [F,pTrust,pArtifact]=qualDiscrimCalls_v2(T,E,homPos,inputParam)
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

aIdx=homPos & T.ACountF+T.ACountR >= T.BCountF+T.BCountR;
bIdx=homPos & T.ACountF+T.ACountR < T.BCountF+T.BCountR;
%%%Create Quality Metrics Table
F=table();
F.TumorPerPassReads=T.ReadDepthPass./T.ReadDepth;
F.normalPerReadPass=T.perReadPass;
F.normalPerReadPass(~isfinite(T.perReadPass))=inputParam.minPerReadPASS;
F.ABfrac=(T.AcountsComb+T.BcountsComb)./T.ReadDepthPass;
F.ABfrac(T.ReadDepthPass==0)=0;
F.normalABfrac=max(T.abFrac,0);
F.minPerStrand=NaN(size(F,1),1);
F.minPerStrand(~homPos)=min([T.ACountF(~homPos)./(T.ACountF(~homPos)+T.ACountR(~homPos)) T.ACountR(~homPos)./(T.ACountF(~homPos)+T.ACountR(~homPos)) T.BCountF(~homPos)./(T.BCountF(~homPos)+T.BCountR(~homPos)) T.BCountR(~homPos)./(T.BCountF(~homPos)+T.BCountR(~homPos))],[],2);
F.minPerStrand(aIdx)=min([T.ACountF(aIdx)./(T.ACountF(aIdx)+T.ACountR(aIdx)) T.ACountR(aIdx)./(T.ACountF(aIdx)+T.ACountR(aIdx))],[],2);
F.minPerStrand(bIdx)=min([T.BCountF(bIdx)./(T.BCountF(bIdx)+T.BCountR(bIdx)) T.BCountR(bIdx)./(T.BCountF(bIdx)+T.BCountR(bIdx))],[],2);
F.minBQ=NaN(size(F,1),1);
F.minBQ(~homPos)=min([T.AmeanBQ(~homPos) T.BmeanBQ(~homPos)],[],2);
F.minBQ(aIdx)=T.AmeanBQ(aIdx);
F.minBQ(bIdx)=T.BmeanBQ(bIdx);
F.minMQ=NaN(size(F,1),1);
F.minMQ(~homPos)=min([T.AmeanMQ(~homPos) T.BmeanMQ(~homPos)],[],2);
F.minMQ(aIdx)=T.AmeanMQ(aIdx);
F.minMQ(bIdx)=T.BmeanMQ(bIdx);
F.maxPMM=NaN(size(F,1),1);
F.maxPMM(~homPos)=max([T.AmeanPMM(~homPos) T.BmeanPMM(~homPos)],[],2);
F.maxPMM(aIdx)=T.AmeanPMM(aIdx);
F.maxPMM(bIdx)=T.BmeanPMM(bIdx);
F.seqEndDist=NaN(size(F,1),1);
F.seqEndDist(~homPos)=min([T.AmeanReadPos(~homPos) inputParam.ReadLength-T.AmeanReadPos(~homPos) T.BmeanReadPos(~homPos) inputParam.ReadLength-T.BmeanReadPos(~homPos)],[],2);   
F.seqEndDist(aIdx)=min([T.AmeanReadPos(aIdx) inputParam.ReadLength-T.AmeanReadPos(aIdx)],[],2);   
F.seqEndDist(bIdx)=min([T.BmeanReadPos(bIdx) inputParam.ReadLength-T.BmeanReadPos(bIdx)],[],2);
F.strandDiff=zeros(size(F,1),1);
F.strandDiff(~homPos)=max(abs(T.ACountF(~homPos)./(T.ACountF(~homPos)+T.ACountR(~homPos))-T.BCountF(~homPos)./(T.BCountF(~homPos)+T.BCountR(~homPos))),0);
F.BQdiff=zeros(size(F,1),1);
F.BQdiff(~homPos)=max(abs(T.AmeanBQ(~homPos)-T.BmeanBQ(~homPos)),0);
F.MQdiff=zeros(size(F,1),1);
F.MQdiff(~homPos)=max(abs(T.AmeanMQ(~homPos)-T.BmeanMQ(~homPos)),0);
F.PMMdiff=zeros(size(F,1),1);
F.PMMdiff(~homPos)=max(abs(T.AmeanPMM(~homPos)-T.BmeanPMM(~homPos)),0);
F.ReadPosDiff=zeros(size(F,1),1);
F.ReadPosDiff(~homPos)=max(abs((T.AmeanReadPos(~homPos)-T.BmeanReadPos(~homPos))),0);
F.posMapQC=-10*log10(T.PosMapQC+1E-6);
F.posMapQC(T.PosMapQC<0)=0;
idx=getPosInRegionSplit([T.Chr T.Pos],[E.Chr E.StartPos E.EndPos+1],inputParam.blockSize);
E.MapQC(E.MapQC<0)=1-1E6;
F.exonMapQC=-10*log10(E.MapQC(idx)+1E-6);
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
    'mean position quality score in exon'};


indelPos=T.A>4 | T.B>4;

%%%find positions that pass strict quality thresholds
goodPos(:,1)=F.TumorPerPassReads>inputParam.minPerReadPASS;
goodPos(:,2)=F.normalPerReadPass>inputParam.minPerReadPASS;
goodPos(:,3)=F.ABfrac>inputParam.minABFrac;
goodPos(:,4)=F.normalABfrac>inputParam.minABFrac;
goodPos(:,5)=F.minPerStrand>inputParam.minPercentStrand;
goodPos(:,6)=F.minBQ>inputParam.minMeanBQ;
goodPos(:,7)=F.minMQ>inputParam.minMeanMQ;
goodPos(:,8)=F.maxPMM<inputParam.maxPMM;
goodPos(:,9)=F.seqEndDist>inputParam.minSeqEndDist;
goodPos(:,10)=F.strandDiff<inputParam.maxStrandDiff;
goodPos(:,11)=F.BQdiff<inputParam.maxBQdiff;
goodPos(:,12)=F.MQdiff<inputParam.maxMQdiff;
goodPos(:,13)=F.PMMdiff<inputParam.maxPMMdiff;
goodPos(:,14)=F.ReadPosDiff<inputParam.maxReadPosDiff;
goodPos(:,15)=F.posMapQC>inputParam.minPosQual;
goodPos(:,16)=F.exonMapQC>inputParam.minExonQual;


%%%find positions that fail minimal quality thresholds
badPos(:,1)=F.TumorPerPassReads<inputParam.perPassReadReject;
badPos(:,2)=F.normalPerReadPass<inputParam.perPassReadReject;
badPos(:,3)=F.ABfrac<inputParam.ABfracReject;
badPos(:,4)=F.normalABfrac<inputParam.ABfracReject;
badPos(:,5)=F.minPerStrand<inputParam.perStrandReject;
badPos(:,6)=F.minBQ<inputParam.meanBQReject;
badPos(:,7)=F.minMQ<inputParam.meanMQReject;
badPos(:,8)=F.maxPMM>inputParam.PMMReject;
badPos(:,9)=F.seqEndDist<inputParam.seqEndDistReject;
badPos(:,10)=F.strandDiff>inputParam.strandDiffReject;
badPos(:,11)=F.BQdiff>inputParam.BQdiffReject;
badPos(:,12)=F.MQdiff>inputParam.MQdiffReject;
badPos(:,13)=F.PMMdiff>inputParam.PMMdiffReject;
badPos(:,14)=F.ReadPosDiff>inputParam.ReadPosDiffReject;
badPos(:,15)=F.posMapQC<inputParam.posQualReject;
badPos(:,16)=F.exonMapQC<inputParam.exonQualReject;

%%%classify SNVs as variants vs artifacts
sampleSNV=F{~indelPos & ~homPos,:};
trainingSNV=[F{sum(goodPos,2)>12 & sum(badPos,2)==0 & ~indelPos & ~homPos,:}; F{sum(badPos,2)>2 & ~indelPos & ~homPos,:}];
groupSNV=[ones(sum(sum(goodPos,2)>12 & sum(badPos,2)==0 & ~indelPos & ~homPos),1); zeros(sum(sum(badPos,2)>2 & ~indelPos & ~homPos),1)];
discrSNV=fitcdiscr(trainingSNV,groupSNV,'DiscrimType', 'pseudoQuadratic');
[~,pArtifact(~indelPos & ~homPos,:)] = predict(discrSNV,sampleSNV);

%%%classify SNVs as trusted vs low quality positions
sampleSNV=F{~indelPos & ~homPos,:};
trainingSNV=[F{sum(goodPos,2)==16 & ~indelPos & ~homPos,:}; F{sum(badPos,2)>0 & ~indelPos & ~homPos,:}];
groupSNV=[ones(sum(sum(goodPos,2)==16 & ~indelPos & ~homPos),1); zeros(sum(sum(badPos,2)>0 & ~indelPos & ~homPos),1)];
discrSNV=fitcdiscr(trainingSNV,groupSNV,'DiscrimType', 'pseudoQuadratic');
[~,pTrust(~indelPos & ~homPos,:)] = predict(discrSNV,sampleSNV);

%%%classify SNVs as variants vs artifacts
sampleSNV=F{~indelPos & homPos,:};
trainingSNV=[F{sum(goodPos,2)>12 & sum(badPos,2)==0 & ~indelPos & homPos,:}; F{sum(badPos,2)>2 & ~indelPos & homPos,:}];
groupSNV=[ones(sum(sum(goodPos,2)>12 & sum(badPos,2)==0 & ~indelPos & homPos),1); zeros(sum(sum(badPos,2)>2 & ~indelPos & homPos),1)];
discrSNV=fitcdiscr(trainingSNV,groupSNV,'DiscrimType', 'pseudoQuadratic');
[~,pArtifact(~indelPos & homPos,:)] = predict(discrSNV,sampleSNV);

%%%classify SNVs as trusted vs low quality positions
sampleSNV=F{~indelPos & homPos,:};
trainingSNV=[F{sum(goodPos,2)==16 & ~indelPos & homPos,:}; F{sum(badPos,2)>0 & ~indelPos & homPos,:}];
groupSNV=[ones(sum(sum(goodPos,2)==16 & ~indelPos & homPos),1); zeros(sum(sum(badPos,2)>0 & ~indelPos & homPos),1)];
discrSNV=fitcdiscr(trainingSNV,groupSNV,'DiscrimType', 'pseudoQuadratic');
[~,pTrust(~indelPos & homPos,:)] = predict(discrSNV,sampleSNV);

%%%classify indels as variants vs artifacts
finitePos=sum(isfinite(F{:,:}),2)==16;
sampleIndel=F{indelPos & finitePos,:};
trainingIndel=[F{sum(goodPos(:,[1:10 12:16]),2)==15 & indelPos & finitePos,:}; F{sum(badPos(:,[1:10 12:16]),2)>0 & indelPos & finitePos,:}];
groupIndel=[ones(sum(sum(goodPos(:,[1:10 12:16]),2)==15 & indelPos & finitePos),1); zeros(sum(sum(badPos(:,[1:10 12:16]),2)>0 & indelPos & finitePos),1)];
discrIndel=fitcdiscr(trainingIndel(:,[1:10 12:16]),groupIndel,'DiscrimType', 'pseudoQuadratic');
[~,pTrust(indelPos & finitePos,:)] = predict(discrIndel,sampleIndel(:,[1:10 12:16]));

%%%classify indels as trusted vs low quality positions
finitePos=sum(isfinite(F{:,:}),2)==16;
sampleIndel=F{indelPos & finitePos,:};
trainingIndel=[F{sum(goodPos(:,[1:10 12:16]),2)>11 & sum(badPos,2)==0 &indelPos & finitePos,:}; F{sum(badPos(:,[1:10 12:16]),2)>2 & indelPos & finitePos,:}];
groupIndel=[ones(sum(sum(goodPos(:,[1:10 12:16]),2)>11 & sum(badPos,2)==0 &indelPos & finitePos),1); zeros(sum(sum(badPos(:,[1:10 12:16]),2)>2 & indelPos & finitePos),1)];
discrIndel=fitcdiscr(trainingIndel(:,[1:10 12:16]),groupIndel,'DiscrimType', 'pseudoQuadratic');
[~,pArtifact(indelPos & finitePos,:)] = predict(discrIndel,sampleIndel(:,[1:10 12:16]));

pTrust(indelPos & ~finitePos,:)=0;
pArtifact(indelPos & ~finitePos,:)=0;

F.goodPosSum=sum(goodPos,2);
F.badPosSum=sum(badPos,2);

F.Properties.VariableDescriptions(17)={'Number of PASS criteria met'};
F.Properties.VariableDescriptions(18)={'Number of REJECT criteria met'};
return;