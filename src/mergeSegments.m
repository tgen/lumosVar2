function mergeSeg=mergeSegments(segs,Ecell,inputParam)

segPosList=[];
for i=1:length(segs)
    segPosList=[segPosList; segs{i}(:,1:2); segs{i}(end,[1 3])];
end

segPosList=sortrows(segPosList,[1 2]);
segPosList=unique(segPosList,'rows');

chrList=unique(segPosList(:,1));
segsMerged=[];
for i=1:length(chrList)
    idx=find(segPosList(:,1)==chrList(i));
    segsMerged=[segsMerged; [ones(length(idx)-1,1)*chrList(i) segPosList(idx(1:end-1),2) segPosList(idx(2:end),2) ones(length(idx)-1,1)]];
end

mergeSeg=removeSegs(segsMerged,Ecell,inputParam,0);