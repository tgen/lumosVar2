function [countsAll]=getCounts(Tcell, inputParam)


posList=[];
for i=1:size(Tcell,2)
    T=Tcell{i};
    posList=[posList; T.Chr T.Pos];
end

posList=unique(posList,'rows');

Amat=NaN(size(posList,1),size(Tcell,2));
Bmat=NaN(size(posList,1),size(Tcell,2));

for i=1:size(Tcell,2)
    T=Tcell{i};
    [lia,locb]=ismember([T.Chr T.Pos],posList,'rows');
    Amat(locb,i)=T.A;
    Bmat(locb,i)=T.B;
    AcountsOrig(locb,i)=T.ACountF+T.ACountR;
    BcountsOrig(locb,i)=T.BCountF+T.BCountR;
    ApopAForig(locb,i)=T.ApopAF;
    BpopAForig(locb,i)=T.BpopAF;
    RDmat(locb,i)=T.ReadDepthPass;
    RefOrig(locb,i)=T.Ref;
    cosmic(locb,i)=T.CosmicCount;
    cnaF(locb,i)=T.cnaF;
    N(locb,i)=T.NumCopies;
    M(locb,i)=T.MinAlCopies;
    Wmat(locb,i)=T.W;
    BmeanBQ(locb,i)=T.BmeanBQ;
    AmeanBQ(locb,i)=T.AmeanBQ;
    mapQC(locb,i)=T.PosMapQC;
end

Ref=zeros(size(Amat));
Ref(min(RefOrig,[],2)==max(RefOrig,[],2),:)=RefOrig(min(RefOrig,[],2)==max(RefOrig,[],2),:);

for i=1:size(Wmat,2)
    Wmat(Wmat(:,i)==0,i)=mode(Wmat(Wmat(:,i)~=0,i));
    BmeanBQ(BmeanBQ(:,i)==0,i)=inputParam.defaultBQ;
    AmeanBQ(AmeanBQ(:,i)==0,i)=inputParam.defaultBQ;
end

A=NaN(size(posList,1),1);
A(min(Amat,[],2)-max(Amat,[],2)==0,:)=Amat(min(Amat,[],2)-max(Amat,[],2)==0,1);
Acounts=zeros(size(Amat));
Acounts(min(Amat,[],2)-max(Amat,[],2)==0,:)=AcountsOrig(min(Amat,[],2)-max(Amat,[],2)==0,:);
ApopAF=zeros(size(Amat));
ApopAF(min(Amat,[],2)-max(Amat,[],2)==0,:)=ApopAForig(min(Amat,[],2)-max(Amat,[],2)==0,:);

B=NaN(size(posList,1),1);
B(min(Bmat,[],2)-max(Bmat,[],2)==0,:)=Bmat(min(Bmat,[],2)-max(Bmat,[],2)==0,1);
Bcounts=zeros(size(Bmat));
Bcounts(min(Bmat,[],2)-max(Bmat,[],2)==0,:)=BcountsOrig(min(Bmat,[],2)-max(Bmat,[],2)==0,:);
BpopAF=zeros(size(Amat));
BpopAF(min(Bmat,[],2)-max(Bmat,[],2)==0,:)=BpopAForig(min(Bmat,[],2)-max(Bmat,[],2)==0,:);

pos=find(isnan(sum([A B Amat Bmat],2)));
for i=1:size(pos)
    %posList(pos(i),:)
    alleles=unique([Amat(pos(i),:) Bmat(pos(i),:)]);
    counts=zeros(size(alleles));
    for j=1:size(alleles,2)
        counts(j)=sum(AcountsOrig(pos(i),alleles(j)==Amat(pos(i),:)))+sum(BcountsOrig(pos(i),alleles(j)==Bmat(pos(i),:)));
    end
    [countSort,idx]=sort(counts,'descend');
    AinsIdx=RefOrig(pos(i),:)<=4 & Amat(pos(i),:)>4;
    BinsIdx=RefOrig(pos(i),:)<=4 & Bmat(pos(i),:)>4;
    BdelIdx=RefOrig(pos(i),:)>4 & Bmat(pos(i),:)<=4;
    AdelIdx=RefOrig(pos(i),:)>4 & Amat(pos(i),:)<=4;
    ArefIdx=Amat(pos(i),:)==RefOrig(pos(i),:);
    BrefIdx=Bmat(pos(i),:)==RefOrig(pos(i),:);
    AaltIdx=Amat(pos(i),:)~=RefOrig(pos(i),:) & Amat(pos(i),:)<=4 & ~AdelIdx;
    BaltIdx=Bmat(pos(i),:)~=RefOrig(pos(i),:) & Bmat(pos(i),:)<=4 & ~BdelIdx;
    insCount=sum([AcountsOrig(pos(i),AinsIdx) BcountsOrig(pos(i),BinsIdx)]);
    delCount=sum([AcountsOrig(pos(i),AdelIdx) BcountsOrig(pos(i),BdelIdx)]);
    refCount=sum([AcountsOrig(pos(i),ArefIdx) BcountsOrig(pos(i),BrefIdx)]);
    altCount=sum([AcountsOrig(pos(i),AaltIdx) BcountsOrig(pos(i),BaltIdx)]);
    [countsSort,tIdx]=sort([refCount altCount insCount delCount],'descend');
    if(tIdx(1)==1)
        refAlleles=unique([Amat(pos(i),ArefIdx) Bmat(pos(i),BrefIdx)]);
        if length(refAlleles)==1;
            A(pos(i),:)=refAlleles;
        elseif tIdx(2)==4 && length(refAlleles(refAlleles>4))==1
            A(pos(i),:)=refAlleles(refAlleles>4);
        elseif tIdx(2)~=4 &&  length(refAlleles(refAlleles<=4))==1
            A(pos(i),:)=refAlleles(refAlleles<=4);
        else
            A(pos(i),:)=-1;
        end
        Acounts(pos(i),:)=AcountsOrig(pos(i),:).*ArefIdx+BcountsOrig(pos(i),:).*BrefIdx;
        ApopAF(pos(i),:)=ApopAForig(pos(i),:).*ArefIdx+BpopAForig(pos(i),:).*BrefIdx;
    else
        if(tIdx(1)==2)
            Aallele=unique([Amat(pos(i),AaltIdx) Bmat(pos(i),BaltIdx)]);
        elseif(tIdx(1)==3)
            Aallele=unique([Amat(pos(i),AinsIdx) Bmat(pos(i),BinsIdx)]);
        else
            Aallele=unique([Amat(pos(i),AdelIdx) Bmat(pos(i),BdelIdx)]);
        end
        if length(Aallele)==1;
            A(pos(i),:)=Aallele;
        else
            lia=ismember(alleles(idx),Aallele);
            A(pos(i),:)=alleles(idx(find(lia,1)));
        end
        aIdx=Amat(pos(i),:)==A(pos(i),:) & ~ArefIdx;
        bIdx=Bmat(pos(i),:)==A(pos(i),:) & ~BrefIdx;
        Acounts(pos(i),:)=AcountsOrig(pos(i),:).*aIdx+BcountsOrig(pos(i),:).*bIdx;
        ApopAF(pos(i),:)=ApopAForig(pos(i),:).*aIdx+BpopAForig(pos(i),:).*bIdx;
    end
    if(tIdx(2)==1)
        refAlleles=unique([Amat(pos(i),ArefIdx) Bmat(pos(i),BrefIdx)]);
        if length(refAlleles)==1;
            B(pos(i),:)=refAlleles;
        elseif tIdx(1)==4 && length(refAlleles(refAlleles>4))==1
            B(pos(i),:)=refAlleles(refAlleles>4);
        elseif tIdx(1)~=4 &&  length(refAlleles(refAlleles<=4))==1
            B(pos(i),:)=refAlleles(refAlleles<=4);
        else
            B(pos(i),:)=-1;
        end
        Bcounts(pos(i),:)=AcountsOrig(pos(i),:).*ArefIdx+BcountsOrig(pos(i),:).*BrefIdx;
        BpopAF(pos(i),:)=ApopAForig(pos(i),:).*ArefIdx+BpopAForig(pos(i),:).*BrefIdx;
    else
        if(tIdx(2)==2)
            Ballele=unique([Amat(pos(i),AaltIdx) Bmat(pos(i),BaltIdx)]);
        elseif(tIdx(2)==3)
            Ballele=unique([Amat(pos(i),AinsIdx) Bmat(pos(i),BinsIdx)]);
        else
            Ballele=unique([Amat(pos(i),AdelIdx) Bmat(pos(i),BdelIdx)]);
        end
        if length(Ballele)==1;
            B(pos(i),:)=Ballele;
        else
            lia=ismember(alleles(idx),Ballele);
            B(pos(i),:)=alleles(idx(find(lia,1)));
        end
        aIdx=Amat(pos(i),:)==B(pos(i),:) & ~ArefIdx;
        bIdx=Bmat(pos(i),:)==B(pos(i),:) & ~BrefIdx;
        Bcounts(pos(i),:)=AcountsOrig(pos(i),:).*aIdx+BcountsOrig(pos(i),:).*bIdx;
        BpopAF(pos(i),:)=ApopAForig(pos(i),:).*aIdx+BpopAForig(pos(i),:).*bIdx;
    end
    if(tIdx(1)==1 && A(pos(i),:)>0)
        Ref(pos(i),:)=A(pos(i),:);
    elseif(tIdx(2)==1 && B(pos(i),:)>0)
        Ref(pos(i),:)=B(pos(i),:);
    elseif(max(tIdx(1:2))==4)
        Ref(pos(i),:)=RefOrig(pos(i),:).*(AdelIdx | BdelIdx);
    else
        Ref(pos(i),:)=RefOrig(pos(i),:).*~(AdelIdx | BdelIdx);
    end
    %i
end

Ref=max(Ref,[],2);
ApopAF=max(ApopAF,[],2);
BpopAF=max(BpopAF,[],2);
cosmic=max(cosmic,[],2);

countsAll=array2table([posList Ref A B ApopAF BpopAF],'VariableNames',{'Chr','Pos','Ref','A','B','ApopAF','BpopAF'});
for i=1:size(Acounts,2)
    countsAll.Acounts(:,i)=Acounts(:,i);
    countsAll.Bcounts(:,i)=Bcounts(:,i);
end
    

            