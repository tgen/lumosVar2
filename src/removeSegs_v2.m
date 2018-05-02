function segsMerged=removeSegs_v2(segsMerged,exonRD,hetData,inputParam,c)


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
    pos2=getPosInRegions(hetData{1}{:,1:2},currRegion);
    for i=1:length(exonRD)
        mat1(:,i)=exonRD{i}(pos==1,4)./exonRD{i}(pos==1,5);
        mat2(:,i)=exonRD{i}(pos==2,4)./exonRD{i}(pos==2,5);
        mat3(:,i)=exonRD{i}(pos==3,4)./exonRD{i}(pos==3,5);
        baf1(:,i)=hetData{i}.BcountsComb(pos2==1)./hetData{i}.ReadDepthPass(pos2==1);
        baf2(:,i)=hetData{i}.BcountsComb(pos2==2)./hetData{i}.ReadDepthPass(pos2==2);
        baf3(:,i)=hetData{i}.BcountsComb(pos2==3)./hetData{i}.ReadDepthPass(pos2==3);
    end
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
segsMerged=removeSegs_v2(segsMerged,exonRD,hetData,inputParam,c);
end