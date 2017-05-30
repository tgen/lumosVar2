function [F,pTrust]=qualDiscrim_v2(T,E,inputParam)
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

%%%Create Quality Metrics Table
F=table();
F.TumorPerPassReads=T.ReadDepthPass./T.ReadDepth;
F.normalPerReadPass=T.perReadPass;
F.normalPerReadPass(~isfinite(T.perReadPass))=inputParam.minPerReadPASS;
F.ABfrac=(T.ACountF+T.ACountR+T.BCountF+T.BCountR)./T.ReadDepthPass;
F.ABfrac(T.ReadDepthPass==0)=0;
F.normalABfrac=max(T.abFrac,0);
F.minPerStrand=min([T.ACountF./(T.ACountF+T.ACountR) T.ACountR./(T.ACountF+T.ACountR) T.BCountF./(T.BCountF+T.BCountR) T.BCountR./(T.BCountF+T.BCountR)],[],2);
F.minBQ=min([T.AmeanBQ T.BmeanBQ],[],2);
F.minMQ=min([T.AmeanMQ T.BmeanMQ],[],2);
F.maxPMM=max([T.AmeanPMM T.BmeanPMM],[],2);
F.seqEndDist=min([T.AmeanReadPos inputParam.ReadLength-T.AmeanReadPos T.BmeanReadPos inputParam.ReadLength-T.BmeanReadPos],[],2);   
F.strandDiff=max(abs(T.ACountF./(T.ACountF+T.ACountR)-T.BCountF./(T.BCountF+T.BCountR)),0);
F.BQdiff=max(abs(T.AmeanBQ-T.BmeanBQ),0);
F.MQdiff=max(abs(T.AmeanMQ-T.BmeanMQ),0);
F.PMMdiff=max(abs(T.AmeanPMM-T.BmeanPMM),0);
F.ReadPosDiff=max(abs((T.AmeanReadPos-T.BmeanReadPos)),0);
F.posMapQC=-10*log10(T.PosMapQC+1E-6);
F.posMapQC(T.PosMapQC<0)=60;
F.posMapQC(isnan(T.PosMapQC))=0;
idx=getPosInRegionSplit([T.Chr T.Pos],[E.Chr E.StartPos E.EndPos],inputParam.blockSize);
F.exonMapQC=-10*log10(E.MapQC(idx)+1E-6);
F.exonMapQC(E.MapQC(idx)<0)=60;
F.exonMapQC(isnan(E.MapQC(idx)))=0;
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

goodPos=F.posMapQC>inputParam.minPosQual & min([T.ApopAF T.BpopAF],[],2)>inputParam.minHetPopFreq & T.BCountF+T.BCountR>inputParam.minBCount;
badPos=F.posMapQC<inputParam.posQualReject & min([T.ApopAF T.BpopAF],[],2)<inputParam.minHetPopFreq & T.BCountF+T.BCountR>inputParam.minBCount;

%%%classify SNVs as variants vs artifacts
sampleSNV=F{~indelPos,:};
trainingSNV=[F{goodPos & ~indelPos,:}; F{badPos & ~indelPos,:}];
groupSNV=[ones(sum(goodPos & ~indelPos),1); zeros(sum(badPos & ~indelPos),1)];
discrSNV=fitcdiscr(trainingSNV,groupSNV,'DiscrimType', 'pseudoQuadratic');
[class(~indelPos),pTrust(~indelPos,:)] = predict(discrSNV,sampleSNV);

%%%classify indels as variants vs artifacts
finitePos=sum(isfinite(F{:,:}),2)==16;
sampleIndel=F{indelPos & finitePos,:};
trainingIndel=[F{goodPos & indelPos & finitePos,:}; F{badPos & indelPos & finitePos,:}];
groupIndel=[ones(sum(goodPos & indelPos & finitePos),1); zeros(sum(badPos & indelPos & finitePos),1)];
discrIndel=fitcdiscr(trainingIndel,groupIndel,'DiscrimType', 'pseudoQuadratic');
[class(indelPos & finitePos),pTrust(indelPos & finitePos,:)] = predict(discrIndel,sampleIndel);

%%%classify indels as trusted vs low quality positions

F.goodPosSum=sum(goodPos,2);
F.badPosSum=sum(badPos,2);

F.Properties.VariableDescriptions(17)={'Number of PASS criteria met'};
F.Properties.VariableDescriptions(18)={'Number of REJECT criteria met'};
return;