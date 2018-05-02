function mergeSeg=mergeSegments(segs,exonRD,Tcell,hetPos,inputParam)

segPosList=[];
for i=1:length(segs)
    segPosList=[segPosList; segs{i}(:,1:2); segs{i}(:,[1 3])];
end

for i=1:length(Tcell)
    hetData{i}=Tcell{i}(hetPos,:);
end

segPosList=sortrows(segPosList,[1 2]);
segPosList=unique(segPosList,'rows');

chrList=unique(segPosList(:,1));
for i=1:length(chrList)
    idx=find(segPosList(:,1)==chrList(i));
    segsMerged=[ones(length(idx)-1,1)*chrList(i) segPosList(idx(1:end-1),2) segPosList(idx(2:end),2)];
    idx=getPosInRegions([exonRD{1}(:,1) mean(exonRD{1}(:,2:3),2)],segsMerged);
    idx2=getPosInRegions(hetData{1}{:,1:2},segsMerged);
    p=nan(size(segsMerged,1)-1,length(exonRD));
    pBaf=nan(size(segsMerged,1)-1,length(exonRD));
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
    segsMergedChr{i}=removeSegs(segsMerged,exonRD,hetData,inputParam,0);
end

mergeSeg=[];
for i=1:length(chrList)
    mergeSeg=[mergeSeg; segsMergedChr{i}];
end
%mergeSeg=removeSegs(segsMerged,Ecell,inputParam,0);
