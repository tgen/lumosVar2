function [F,pTrust,pArtifact]=qualDiscrimCalls(T,E,homPos,inputParam)
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
%    homPos - logical vector of positions that are called homzygous
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
% See also: LumosVarMain, getCounts, jointSNV

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

%numVar=gather(height(T));

%%%set quality metrics of positions with no B allele to NaN
T.BmeanBQ(T.BCountF+T.BCountR==0)=NaN;
T.BmeanMQ(T.BCountF+T.BCountR==0)=NaN;
T.BmeanPMM(T.BCountF+T.BCountR==0)=NaN;
T.BmeanReadPos(T.BCountF+T.BCountR==0)=NaN;

%%%determine whether homozygous positions are A or B
aIdx=homPos & T.ACountF+T.ACountR >= T.BCountF+T.BCountR;
bIdx=homPos & T.ACountF+T.ACountR < T.BCountF+T.BCountR;

readLength=max([T.AmeanReadPos; T.BmeanReadPos]);
inputParam.minSeqEndDist=readLength*(inputParam.minSeqEndDist./inputParam.ReadLength);
inputParam.ReadLength=readLength;

%%%Create Quality Metrics Table
F=T(:,2:3);
F.TumorPerPassReads=(T.ReadDepthPass+1)./(T.ReadDepth+1);
F.normalPerReadPass=T.perReadPass;
F.normalPerReadPass(~isfinite(T.perReadPass))=inputParam.minPerReadPASS;
F.ABfrac=(T.AcountsComb+T.BcountsComb+1)./(T.ReadDepthPass+1);
F.normalABfrac=T.abFrac;
F.normalABfrac(~isfinite(T.abFrac))=inputParam.minABFrac;
F.minPerStrand=min(min([T.ACountF./(T.ACountF+T.ACountR) T.ACountR./(T.ACountF+T.ACountR) T.BCountF./(T.BCountF+T.BCountR) T.BCountR./(T.BCountF+T.BCountR)],[],2),0.5);
%T.minPerStrand(~homPos)=min(min([T.ACountF(~homPos)./(T.ACountF(~homPos)+T.ACountR(~homPos)) T.ACountR(~homPos)./(T.ACountF(~homPos)+T.ACountR(~homPos)) T.BCountF(~homPos)./(T.BCountF(~homPos)+T.BCountR(~homPos)) T.BCountR(~homPos)./(T.BCountF(~homPos)+T.BCountR(~homPos))],[],2),0.5);
F.minPerStrand(aIdx)=min(min([T.ACountF(aIdx)./(T.ACountF(aIdx)+T.ACountR(aIdx)) T.ACountR(aIdx)./(T.ACountF(aIdx)+T.ACountR(aIdx))],[],2),0.5);
F.minPerStrand(bIdx)=min(min([T.BCountF(bIdx)./(T.BCountF(bIdx)+T.BCountR(bIdx)) T.BCountR(bIdx)./(T.BCountF(bIdx)+T.BCountR(bIdx))],[],2),0.5);
F.minBQ=min([T.AmeanBQ T.BmeanBQ],[],2);
F.minBQ(aIdx)=T.AmeanBQ(aIdx);
F.minBQ(bIdx)=T.BmeanBQ(bIdx);
F.minMQ=min([T.AmeanMQ T.BmeanMQ],[],2);
F.minMQ(aIdx)=T.AmeanMQ(aIdx);
F.minMQ(bIdx)=T.BmeanMQ(bIdx);
F.maxPMM=max([T.AmeanPMM T.BmeanPMM],[],2);
F.maxPMM(aIdx)=T.AmeanPMM(aIdx);
F.maxPMM(bIdx)=T.BmeanPMM(bIdx);
F.seqEndDist=min([T.AmeanReadPos inputParam.ReadLength-T.AmeanReadPos T.BmeanReadPos inputParam.ReadLength-T.BmeanReadPos],[],2);   
F.seqEndDist(aIdx)=min([T.AmeanReadPos(aIdx) inputParam.ReadLength-T.AmeanReadPos(aIdx)],[],2);   
F.seqEndDist(bIdx)=min([T.BmeanReadPos(bIdx) inputParam.ReadLength-T.BmeanReadPos(bIdx)],[],2);
F.strandDiff=max(abs(T.ACountF./(T.ACountF+T.ACountR)-T.BCountF./(T.BCountF+T.BCountR)),0);
F.strandDiff(homPos)=0;
F.BQdiff=max(abs(T.AmeanBQ-T.BmeanBQ),0);
F.BQdiff(homPos)=0;
F.MQdiff=max(abs(T.AmeanMQ-T.BmeanMQ),0);
F.MQdiff(homPos)=0;
F.PMMdiff=max(abs(T.AmeanPMM-T.BmeanPMM),0);
F.PMMdiff(homPos)=0;
F.ReadPosDiff=max(abs((T.AmeanReadPos-T.BmeanReadPos)),0);
F.ReadPosDiff(homPos)=0;
F.posMapQC=-10*log10(T.PosMapQC+1E-6);
F.posMapQC(T.PosMapQC<0)=0;
F.posMapQC(~isfinite(T.PosMapQC))=inputParam.minPosQual;
%idx=getPosInRegionSplit([T.Chr T.Pos],[E.Chr E.StartPos E.EndPos+1],inputParam.blockSize);
F.exonMapQC=F.posMapQC;
% E.MapQC(E.MapQC<0)=1-1E6;
% F.exonMapQC=zeros(numVar,1);
% F.exonMapQC(~isnan(idx))=-10*log10(E.MapQC(idx(~isnan(idx)))+1E-6);
F.Properties.VariableDescriptions={'Chr','Pos','Percent of Reads in Tumor Passing Quality Thresh', ...
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
F.goodPos=F.TumorPerPassReads>inputParam.minPerReadPASS;
F.goodPos=[F.goodPos F.normalPerReadPass>inputParam.minPerReadPASS];
F.goodPos=[F.goodPos F.ABfrac>inputParam.minABFrac];
F.goodPos=[F.goodPos F.normalABfrac>inputParam.minABFrac];
F.goodPos=[F.goodPos F.minPerStrand>inputParam.minPercentStrand];
F.goodPos=[F.goodPos F.minBQ>inputParam.minMeanBQ];
F.goodPos=[F.goodPos F.minMQ>inputParam.minMeanMQ];
F.goodPos=[F.goodPos F.maxPMM<inputParam.maxPMM];
F.goodPos=[F.goodPos F.seqEndDist>inputParam.minSeqEndDist];
F.goodPos=[F.goodPos F.strandDiff<inputParam.maxStrandDiff];
F.goodPos=[F.goodPos F.BQdiff<inputParam.maxBQdiff];
F.goodPos=[F.goodPos F.MQdiff<inputParam.maxMQdiff];
F.goodPos=[F.goodPos F.PMMdiff<inputParam.maxPMMdiff];
F.goodPos=[F.goodPos F.ReadPosDiff<inputParam.maxReadPosDiff];
F.goodPos=[F.goodPos F.posMapQC>inputParam.minPosQual];
F.goodPos=[F.goodPos F.exonMapQC>inputParam.minExonQual];


%%%find positions that fail minimal quality thresholds
F.badPos=F.TumorPerPassReads<inputParam.perPassReadReject;
F.badPos=[F.badPos F.normalPerReadPass<inputParam.perPassReadReject];
F.badPos=[F.badPos F.ABfrac<inputParam.ABfracReject];
F.badPos=[F.badPos F.normalABfrac<inputParam.ABfracReject];
F.badPos=[F.badPos F.minPerStrand<inputParam.perStrandReject];
F.badPos=[F.badPos F.minBQ<inputParam.meanBQReject];
F.badPos=[F.badPos F.minMQ<inputParam.meanMQReject];
F.badPos=[F.badPos F.maxPMM>inputParam.PMMReject];
F.badPos=[F.badPos F.seqEndDist<inputParam.seqEndDistReject];
F.badPos=[F.badPos F.strandDiff>inputParam.strandDiffReject];
F.badPos=[F.badPos F.BQdiff>inputParam.BQdiffReject];
F.badPos=[F.badPos F.MQdiff>inputParam.MQdiffReject];
F.badPos=[F.badPos F.PMMdiff>inputParam.PMMdiffReject];
F.badPos=[F.badPos F.ReadPosDiff>inputParam.ReadPosDiffReject];
F.badPos=[F.badPos F.posMapQC<inputParam.posQualReject];
F.badPos=[F.badPos F.exonMapQC<inputParam.exonQualReject];

%%%classify SNVs as variants vs artifacts
F.trainSNVvar=sum(F.goodPos,2)>12 & sum(F.badPos,2)==0 & ~indelPos & ~homPos;
F.trainSNVart=sum(F.badPos,2)>2 & ~indelPos & ~homPos;
F.trainSNV1=F.trainSNVart | F.trainSNVvar;
trainSNV1=[F(F.trainSNV1,2:18) F(F.trainSNV1,"trainSNVvar")];
%groupSNV1=F.trainSNVvar(F.trainSNV1);
sampleSNV=F{:,2:18};
% trainingSNV=[F{sum(F.goodPos,2)>12 & sum(F.badPos,2)==0 & ~indelPos & ~homPos,2:18}; F{sum(F.badPos,2)>2 & ~indelPos & ~homPos,2:18}];
% groupSNV=[true(gather(sum(sum(goodPos,2)>12 & sum(badPos,2)==0 & ~indelPos & ~homPos)),1); false(gather(sum(sum(badPos,2)>2 & ~indelPos & ~homPos)),1)];
discrSNV=fitcdiscr(trainSNV1,'trainSNVvar','DiscrimType', 'pseudoQuadratic');
[~,pArtifact] = predict(discrSNV,sampleSNV);

F.trainSNVpass=sum(F.goodPos,2)==16 & ~indelPos & ~homPos;
F.trainSNVlqc=sum(F.badPos,2)>0 & ~indelPos & ~homPos
F.trainSNV2=F.trainSNVpass | F.trainSNVlqc;
trainSNV2=[F(F.trainSNV2,2:18) F(F.trainSNV2,"trainSNVvar")];
discrSNV2=fitcdiscr(trainSNV2,'trainSNVvar','DiscrimType', 'pseudoQuadratic');
[~,pTrust] = predict(discrSNV2,sampleSNV);

return;
