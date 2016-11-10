function segsMerged=removeSegs(segsMerged,Ecell,inputParam,c)

idx=getPosInRegions([Ecell{1}.Chr mean([Ecell{1}.StartPos Ecell{1}.EndPos],2)],segsMerged);
counts=hist(idx,[1:length(segsMerged)]);
counts(segsMerged(:,4)==0)=Inf;
[~,sortIdx]=sort(counts);

if(sortIdx(1)==1)
    h1=1;
elseif(segsMerged(sortIdx(1)-1,1)==segsMerged(sortIdx(1),1))
   mat1=zeros(sum(idx==sortIdx(1)-1),length(Ecell));
   mat2=zeros(sum(idx==sortIdx(1)),length(Ecell));
   for i=1:length(Ecell)
       mat1(:,i)=Ecell{i}.TumorRD(idx==sortIdx(1)-1)./Ecell{i}.NormalRD(idx==sortIdx(1)-1);
       mat2(:,i)=Ecell{i}.TumorRD(idx==sortIdx(1))./Ecell{i}.NormalRD(idx==sortIdx(1));
   end
   if(size(mat1,1)>1)
       [h1,p1]=ttest2(mat1,mat2,inputParam.cnaAlpha);
   elseif(size(mat1,1)==1)
       [h1,p1]=ttest2(mat2,mat1,inputParam.cnaAlpha);
   else
       h1=0;
       p1=1;
   end
else
   h1=1;
end

if(sortIdx(1)==length(counts))
    h2=1;
elseif(segsMerged(sortIdx(1)+1,1)==segsMerged(sortIdx(1),1))
   mat1=zeros(sum(idx==sortIdx(1)+1),length(Ecell));
   mat2=zeros(sum(idx==sortIdx(1)),length(Ecell));
   for i=1:length(Ecell)
       mat1(:,i)=Ecell{i}.TumorRD(idx==sortIdx(1)+1)./Ecell{i}.NormalRD(idx==sortIdx(1)+1);
       mat2(:,i)=Ecell{i}.TumorRD(idx==sortIdx(1))./Ecell{i}.NormalRD(idx==sortIdx(1));
   end
   if(size(mat1,1)>1)
       [h2,p2]=ttest2(mat1,mat2,inputParam.cnaAlpha);
   elseif(size(mat1,1)==1)
       %size(mat1)
       %size(mat2)
       [h2,p2]=ttest2(mat2,mat1,inputParam.cnaAlpha);
   else
       h2=0;
       p2=1;
   end
else
    h2=1;
end

if(sum(h1)==0 && sum(h2)==0)
    if(min(p1)>min(p2))
        newBound=segsMerged(sortIdx(1),3);
    else
        newBound=segsMerged(sortIdx(1),2);
    end
elseif(sum(h1)==0)
    newBound=segsMerged(sortIdx(1),3);
elseif(sum(h2)==0)
    newBound=segsMerged(sortIdx(1),2);
else
    newBound=NaN;
end

if(isnan(newBound))
    segsMerged(sortIdx(1),4)=0;
else
    if(sortIdx(1)>1)
        if(segsMerged(sortIdx(1)-1,1)==segsMerged(sortIdx(1),1))
            segsMerged(sortIdx(1)-1,3)=newBound;
        end
    end
    if(sortIdx(1)<length(counts))
        if(segsMerged(sortIdx(1)+1,1)==segsMerged(sortIdx(1),1))   
            segsMerged(sortIdx(1)+1,2)=newBound;
        end
    end
    segsMerged=segsMerged([1:sortIdx(1)-1 sortIdx(1)+1:end],:);
end

if(sum(segsMerged(:,4))==0)
    return;
else
    sum(segsMerged(:,4));
    c=c+1;
    segsMerged=removeSegs(segsMerged,Ecell,inputParam,c);
end
    