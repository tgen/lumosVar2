function fStart=makeFstart(fOld)

m=1;
distMat=ones(size(fOld,1),1000);
fStart=[];
while m./max(distMat(:))>0.2
    fRand=100*rand(size(fOld,1),1000);
    distMat=squareform(pdist([fOld fStart fRand]','cityblock'));
    minList=min(distMat(:,1:size(fOld,2)+size(fStart,2)),[],2);
    [m,idx]=max(minList(size(fOld,2)+size(fStart,2)+1:end));
    fStart=[fStart fRand(:,idx)];
end