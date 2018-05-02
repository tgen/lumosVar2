function [segs,bafSegs]=segmentData(exonRD,hetData,cnaAlpha)
%segmentData - parellizes segmentation by chromosome
%uses circular binary segmentation from bioinformatics toolbox
%
% Syntax: segs=segmentData(exonRD,cnaAlpha)
%
% Inputs:
%   exonRD: matrix of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   cnaAlpha: significance cutoff for segmentation
%   
% Outputs:
%   segs: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%
% Other m-files required: none
% Other requirements: bioinformatics toolbox
% Subfunctions: none
% MAT-files required: cghcbshybridnu.mat
%
% See also: TumorOnlyWrapper, cghcbs

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------
numChr=length(unique(exonRD(:,1)));

parfor i=1:numChr
    clustsegs{i}=cghcbs([exonRD(exonRD(:,1)==i,1) exonRD(exonRD(:,1)==i,2) log((exonRD(exonRD(:,1)==i,4)+1)./(exonRD(exonRD(:,1)==i,5)+1))],'Alpha',cnaAlpha);
    bafclustsegs{i}=cghcbs([hetData.Chr(hetData.Chr==i) hetData.Pos(hetData.Chr==i)  hetData.BcountsComb(hetData.Chr==i)./hetData.ReadDepthPass(hetData.Chr==i)],'Alpha',cnaAlpha);
end
    
segs=[];
bafSegs=[];
for i=1:numChr
    segs=[segs; [clustsegs{i}.SegmentData(1).Chromosome*ones(size(clustsegs{i}.SegmentData(1).Start)) clustsegs{i}.SegmentData(1).Start clustsegs{i}.SegmentData(1).End clustsegs{i}.SegmentData(1).Mean]];
    if (~isempty(bafclustsegs{i}))
        bafSegs=[bafSegs; [bafclustsegs{i}.SegmentData(1).Chromosome*ones(size(bafclustsegs{i}.SegmentData(1).Start)) bafclustsegs{i}.SegmentData(1).Start bafclustsegs{i}.SegmentData(1).End bafclustsegs{i}.SegmentData(1).Mean]];
    end
end
return;
