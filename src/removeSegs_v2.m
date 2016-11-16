function segsMerged=removeSegs_v2(segsMerged,exonRD,inputParam,c)


if(sum(segsMerged(:,4)>inputParam.cnaAlpha)==0)
    return
else
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
    pos=getPosInRegions([exonRD{1}(:,1) mean(exonRD{1}(:,2:3),2)],currRegion);
    for i=1:length(exonRD)
        mat1(:,i)=exonRD{i}(pos==1,4)./exonRD{i}(pos==1,5);
        mat2(:,i)=exonRD{i}(pos==2,4)./exonRD{i}(pos==2,5);
        mat3(:,i)=exonRD{i}(pos==3,4)./exonRD{i}(pos==3,5);
    end
    if idx>1
        if size(mat1,1)>1
            [~,p]=ttest2(mat1,mat2,inputParam.cnaAlpha);
            segsMerged(idx-1,4)=min(min(p),1);
        elseif size(mat2,1)>1
            [~,p]=ttest2(mat2,mat1,inputParam.cnaAlpha);
            segsMerged(idx-1,4)=min(min(p),1);
        else
            segsMerged(idx-1,4)=1;
        end
    end
    if idx<size(segsMerged,1)
        if size(mat2,1)>1
            [~,p]=ttest2(mat2,mat3,inputParam.cnaAlpha);
            segsMerged(idx,4)=min(min(p),1);
        elseif size(mat3,1)>1
            [~,p]=ttest2(mat3,mat2,inputParam.cnaAlpha);
            segsMerged(idx,4)=min(min(p),1);
        else
            segsMerged(idx,4)=1;
        end
    end
    c=c+1;
    segsMerged=removeSegs_v2(segsMerged,exonRD,inputParam,c);
end