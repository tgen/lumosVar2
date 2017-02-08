function fStart=makeFstart(fOld,tIdx,inputParam)

m=1;
distMat=ones(size(fOld,1),1000);
fStart=[];
while m./max(distMat(:))>inputParam.minStartDiff
    %fRand=100*rand(size(fOld,1),1000);
    fRand=100*betarnd(ones(length(tIdx),1000)*inputParam.alphaF,((inputParam.alphaF-1)./inputParam.priorF(tIdx)-inputParam.alphaF+2)*ones(1,1000),size(fOld,1),1000);
    distMat=squareform(pdist([fOld fStart fRand]','cityblock'));
    minList=min(distMat(:,1:size(fOld,2)+size(fStart,2)),[],2);
    [m,idx]=max(minList(size(fOld,2)+size(fStart,2)+1:end));
    fStart=[fStart fRand(:,idx)];
end