function [segs,bafSegs]=segmentData(exonRD,hetData,cnaAlpha)
%segmentData - parellizes segmentation by chromosome
%uses circular binary segmentation from bioinformatics toolbox
%
% Syntax: [segs,bafSegs=segmentData(exonRD,hetData,cnaAlpha)
%
% Inputs:
%   exonRD: cell array of matri of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   hetData: table of heterozygous positions with columns: {'Chr', 'Pos', 
%       'BcountsComb','ReadDepthPass'}
%   cnaAlpha: significance cutoff for segmentation
%   
% Outputs:
%   segs: matrix of read depth segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%   bafSegs: matrix of b-allele frequency segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'mean B allele Fraction'
%
% Other m-files required: none
% Other requirements: bioinformatics toolbox
% Subfunctions: none
% MAT-files required: cghcbshybridnu.mat
%
% See also: LumosVarMain, cghcbs

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 9-May-2018

%------------- BEGIN CODE --------------

%%%run segmentation by Chr
numChr=length(unique(exonRD(:,1)));
parfor i=1:numChr
    clustsegs{i}=cghcbs([exonRD(exonRD(:,1)==i,1) exonRD(exonRD(:,1)==i,2) log((exonRD(exonRD(:,1)==i,4)+1)./(exonRD(exonRD(:,1)==i,5)+1))],'Alpha',cnaAlpha);
    bafclustsegs{i}=cghcbs([hetData.Chr(hetData.Chr==i) hetData.Pos(hetData.Chr==i)  hetData.BcountsComb(hetData.Chr==i)./hetData.ReadDepthPass(hetData.Chr==i)],'Alpha',cnaAlpha);
end

%%%combine segments
segs=[];
bafSegs=[];
for i=1:numChr
    segs=[segs; [clustsegs{i}.SegmentData(1).Chromosome*ones(size(clustsegs{i}.SegmentData(1).Start)) clustsegs{i}.SegmentData(1).Start clustsegs{i}.SegmentData(1).End clustsegs{i}.SegmentData(1).Mean]];
    if (~isempty(bafclustsegs{i}))
        bafSegs=[bafSegs; [bafclustsegs{i}.SegmentData(1).Chromosome*ones(size(bafclustsegs{i}.SegmentData(1).Start)) bafclustsegs{i}.SegmentData(1).Start bafclustsegs{i}.SegmentData(1).End bafclustsegs{i}.SegmentData(1).Mean]];
    end
end
return;
