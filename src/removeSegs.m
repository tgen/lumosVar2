function segsMerged=removeSegs(segsMerged,exonRD,hetData,inputParam,c)
%removeSegs - recursively removes non-significant segmentation boundaries
%
% Syntax: segsMerged=removeSegs(segsMerged,exonRD,hetData,inputParam,c)
%
% Inputs:
%   segsMerged: matrix with  the following columns
%       1-'Chr',2-'StartPos',3-'EndPos',4-'boundary geomean p-value'
%   exonRD: cell array of matrices of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   hetData - cell array of tables with length equal to the number of bams,
%       where each row represents a heterozygous variant, and each table 
%       must have the following columns:  {'Chr','Pos', 'ReadDepthPass',
%       'BCountsComb'}
%   inputParam - structure with the folling fields: CNAalpha
%   c - count of recursions
%   
% Outputs:
%   segsMerged - matrix of segmentation boundaries with columns: 1-'Chr',
%           2-'StartPos',3-'EndPos',4-'boundary geomean p-value'
%       
%
% Other m-files required: none
% Other requirements: none
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain, mergeSegments

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 9-May-2018

%------------- BEGIN CODE --------------

%%% exit if all p-values are significant
if(sum(segsMerged(:,4)>inputParam.cnaAlpha)==0)
    return  
else
    %%%remove leaset significant boundary
    [~,idx]=max(segsMerged(:,4));
    segsMerged=[segsMerged(1:idx-1,:); [segsMerged(idx,1:2) segsMerged(idx+1,3) NaN]; segsMerged(idx+2:end,:)];
    currRegion=zeros(3,3);
    if(idx>1)
        currRegion(1,:)=segsMerged(idx-1,1:3);
    end
    currRegion(2,:)=segsMerged(idx,1:3);
    if(idx<size(segsMerged,1))
        currRegion(3,:)=segsMerged(idx+1,1:3);
    end
    %%%get read depth ratio of and b-allele freq of newly merged segment
    %%%and neighbors
    pos=getPosInRegions([exonRD{1}(:,1) mean(exonRD{1}(:,2:3),2)],currRegion);
    pos2=getPosInRegions(hetData{1}{:,1:2},currRegion);
    for i=1:length(exonRD)
        mat1(:,i)=exonRD{i}(pos==1,4)./exonRD{i}(pos==1,5);
        mat2(:,i)=exonRD{i}(pos==2,4)./exonRD{i}(pos==2,5);
        mat3(:,i)=exonRD{i}(pos==3,4)./exonRD{i}(pos==3,5);
        baf1(:,i)=hetData{i}.BcountsComb(pos2==1)./hetData{i}.ReadDepthPass(pos2==1);
        baf2(:,i)=hetData{i}.BcountsComb(pos2==2)./hetData{i}.ReadDepthPass(pos2==2);
        baf3(:,i)=hetData{i}.BcountsComb(pos2==3)./hetData{i}.ReadDepthPass(pos2==3);
    end
    %%%find significance of newly merged segment with neighbors
    p=ones(1,size(mat1,2));
    pBaf=ones(1,size(mat1,2));
    if idx>1
        if size(mat1,1)>1
            [~,p]=ttest2(mat1,mat2,inputParam.cnaAlpha);
        elseif size(mat2,1)>1
            [~,p]=ttest2(mat2,mat1,inputParam.cnaAlpha);
        end
        if size(baf1,1)>1
            [~,pBaf]=ttest2(baf1,baf2,inputParam.cnaAlpha);
        elseif size(baf2,1)>1
            [~,pBaf]=ttest2(baf2,baf1,inputParam.cnaAlpha);
        end
        p=min(p,1);
        pBaf=min(pBaf,1);
        segsMerged(idx-1,4)=min(geomean([p pBaf],2),1);
    end
    if idx<size(segsMerged,1)
        if size(mat2,1)>1
            [~,p]=ttest2(mat2,mat3,inputParam.cnaAlpha);
        elseif size(mat3,1)>1
            [~,p]=ttest2(mat3,mat2,inputParam.cnaAlpha);
        end
        if size(baf2,1)>1
            [~,pBaf]=ttest2(baf2,baf3,inputParam.cnaAlpha);
        elseif size(baf3,1)>1
            [~,pBaf]=ttest2(baf3,baf2,inputParam.cnaAlpha);
        end
        p=min(p,1);
        pBaf=min(pBaf,1);
        segsMerged(idx,4)=min(geomean([p pBaf],2),1);
    end
c=c+1;
%%%try removing more segments
segsMerged=removeSegs(segsMerged,exonRD,hetData,inputParam,c);
end