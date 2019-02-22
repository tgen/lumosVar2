function mergeSeg=mergeSegments(segs,exonRD,hetData,inputParam)
%segmentData - takes segments from different samples and data types and
%merges them
%
% Syntax: mergeSeg=mergeSegments(segs,exonRD,Tcell,hetPos,inputParam)
%
% Inputs:
%   segs: cell array of matrices with  the following columns
%       1-'Chr',2-'StartPos',3-'EndPos',4-'mean of segment'
%   exonRD: cell array of matri of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   Tcell - cell array of tables with length equal to the number of bams,
%       each table must have the following columns:  {'Chr','Pos', 'ReadDepthPass',
%       'BCountsComb'}
%   hetPos - vector of positions in of heterozygous positions
%   inputParam - structure with the folling fields: CNAalpha
%   
% Outputs:
%   mergeSeg - matrix of segmentation boundaries with columns: 1-'Chr',
%           2-'StartPos',3-'EndPos',4-'boundary geomean p-value'
%       
%
% Other m-files required: removeSegs
% Other requirements: none
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain, segmentData

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 9-May-2018

%------------- BEGIN CODE --------------

%%%combine segmentation boundaries
segPosList=[];
for i=1:length(segs)
    segPosList=[segPosList; segs{i}(:,1:2); segs{i}(:,[1 3])];
end
segPosList=sortrows(segPosList,[1 2]);
segPosList=unique(segPosList,'rows');

%%%examine segment data for each chromosome
chrList=unique(segPosList(:,1));
for i=1:length(chrList)
    idx=find(segPosList(:,1)==chrList(i));
    segsMerged=[ones(length(idx)-1,1)*chrList(i) segPosList(idx(1:end-1),2) segPosList(idx(2:end),2)];
    idx=getPosInRegions([exonRD{1}(:,1) mean(exonRD{1}(:,2:3),2)],segsMerged);
    idx2=getPosInRegions(hetData{1}{:,1:2},segsMerged);
    p=nan(size(segsMerged,1)-1,length(exonRD));
    pBaf=nan(size(segsMerged,1)-1,length(exonRD));
    %%%use two sample ttest to find significane of difference in read depth
    %%%ratios and b-allele fractions between neighboring segments
    for j=1:size(segsMerged,1)-1
        mat1=nan(sum(idx==j),length(exonRD));
        mat2=nan(sum(idx==j+1),length(exonRD));
        baf1=nan(sum(idx2==j),length(exonRD));
        baf2=nan(sum(idx2==j+1),length(exonRD));
        for k=1:length(exonRD)
            mat1(:,k)=exonRD{k}(idx==j,4)./exonRD{k}(idx==j,5);
            mat2(:,k)=exonRD{k}(idx==j+1,4)./exonRD{k}(idx==j+1,5);
            baf1(:,k)=hetData{k}.BcountsComb(idx2==j)./hetData{k}.ReadDepthPass(idx2==j);
            baf2(:,k)=hetData{k}.BcountsComb(idx2==j+1)./hetData{k}.ReadDepthPass(idx2==j+1);
        end
        if size(mat1,1)>1
            [~,p(j,:)]=ttest2(mat1,mat2,inputParam.cnaAlpha);
        elseif size(mat2,1)>1
            [~,p(j,:)]=ttest2(mat2,mat1,inputParam.cnaAlpha);
        else
            p(j,:)=1;
        end
        if size(baf1,1)>1
            [~,pBaf(j,:)]=ttest2(baf1,baf2,inputParam.cnaAlpha);
        elseif size(baf2,1)>1
            [~,pBaf(j,:)]=ttest2(baf2,baf1,inputParam.cnaAlpha);
        else
            pBaf(j,:)=1;
        end
    end
    segsMerged(:,4)=[min(geomean([p pBaf],2),1); NaN];
    %%%recursively remove non-significant segment boundaries
    segsMergedChr{i}=removeSegs(segsMerged,exonRD,hetData,inputParam,0);
end

%%%combine chromosome segments
mergeSeg=[];
for i=1:length(chrList)
    mergeSeg=[mergeSeg; segsMergedChr{i}];
end