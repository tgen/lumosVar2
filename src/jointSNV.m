function [postComb, pDataSum, pDataComb,cloneId,prior,countsAll]=jointSNV(Tcell, exonRD, fIn, W, inputParam)

f=zeros(length(Tcell),inputParam.numClones);
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
f(tIdx,1:end)=reshape(fIn,[],inputParam.numClones);

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
    alleles=unique([Amat(pos(i),:) Bmat(pos(i),:)]);
    counts=zeros(size(alleles));
    for j=1:size(alleles,2)
        counts(j)=sum(AcountsOrig(pos(i),alleles(j)==Amat(pos(i),:)))+sum(BcountsOrig(pos(i),alleles(j)==Bmat(pos(i),:)));
    end
    [countSort,idx]=sort(counts,'descend');
    if(max(alleles(idx<=2))>4)
        AinsIdx=Amat(pos(i),:)~=RefOrig(pos(i),:) & Amat(pos(i),:)>4;
        BinsIdx=Bmat(pos(i),:)~=RefOrig(pos(i),:) & Bmat(pos(i),:)>4;
        BdelIdx=Amat(pos(i),:)==RefOrig(pos(i),:) & Amat(pos(i),:)>4;
        AdelIdx=(Bmat(pos(i),:)==RefOrig(pos(i),:) & Bmat(pos(i),:)>4) | (RefOrig(pos(i),:)>4 & BcountsOrig(pos(i),:)==0);
        ArefIdx=Amat(pos(i),:)==RefOrig(pos(i),:);
        BrefIdx=Bmat(pos(i),:)==RefOrig(pos(i),:);
        AaltIdx=Amat(pos(i),:)~=RefOrig(pos(i),:) & Amat(pos(i),:)<=4 & ~AdelIdx;
        BaltIdx=Bmat(pos(i),:)~=RefOrig(pos(i),:) & Bmat(pos(i),:)<=4 & ~BdelIdx;
        insCount=sum([AcountsOrig(pos(i),AinsIdx) BcountsOrig(pos(i),BinsIdx)]);
        delCount=sum([AcountsOrig(pos(i),AdelIdx) BcountsOrig(pos(i),BdelIdx)]);
        refCount=sum([AcountsOrig(pos(i),ArefIdx) BcountsOrig(pos(i),BrefIdx)]);
        altCount=sum([AcountsOrig(pos(i),AaltIdx) BcountsOrig(pos(i),BaltIdx)]);
        [~,tIdx]=sort([refCount altCount insCount delCount],'descend');
        if(tIdx(1)==4)
            if length(unique([Amat(pos(i),AdelIdx) Bmat(pos(i),BdelIdx)]))>1
                A(pos(i),:)=-1;
            else
                A(pos(i),:)=unique([Amat(pos(i),AdelIdx) Bmat(pos(i),BdelIdx)]);
            end
            Acounts(pos(i),:)=AcountsOrig(pos(i),:).*AdelIdx+BcountsOrig(pos(i),:).*BdelIdx;
            ApopAF(pos(i),:)=ApopAForig(pos(i),:).*AdelIdx+BpopAForig(pos(i),:).*BdelIdx;
        elseif(tIdx(1)==3)
            if length(unique([Amat(pos(i),AinsIdx) Bmat(pos(i),BinsIdx)]))>1
                A(pos(i),:)=-1;
            else
                A(pos(i),:)=unique([Amat(pos(i),AinsIdx) Bmat(pos(i),BinsIdx)]);
            end
            Acounts(pos(i),:)=AcountsOrig(pos(i),:).*AinsIdx+BcountsOrig(pos(i),:).*BinsIdx;
            ApopAF(pos(i),:)=ApopAForig(pos(i),:).*AinsIdx+BpopAForig(pos(i),:).*BinsIdx;
        elseif(tIdx(1)==1)
            if length(unique([Amat(pos(i),ArefIdx) Bmat(pos(i),BrefIdx)]))>1
                A(pos(i),:)=-1;
            else
                A(pos(i),:)=unique([Amat(pos(i),ArefIdx) Bmat(pos(i),BrefIdx)]);
            end
            Acounts(pos(i),:)=AcountsOrig(pos(i),:).*ArefIdx+BcountsOrig(pos(i),:).*BrefIdx;
            ApopAF(pos(i),:)=ApopAForig(pos(i),:).*ArefIdx+BpopAForig(pos(i),:).*BrefIdx;
        else
            if length(unique([Amat(pos(i),AaltIdx) Bmat(pos(i),BaltIdx)]))>1
                A(pos(i),:)=-1;
            else
                A(pos(i),:)=unique([Amat(pos(i),AaltIdx) Bmat(pos(i),BaltIdx)]);
            end
            Acounts(pos(i),:)=AcountsOrig(pos(i),:).*AaltIdx+BcountsOrig(pos(i),:).*BaltIdx;
            ApopAF(pos(i),:)=ApopAForig(pos(i),:).*AaltIdx+BpopAForig(pos(i),:).*BaltIdx;
        end
        if(tIdx(2)==4)
            if length(unique([Amat(pos(i),AdelIdx) Bmat(pos(i),BdelIdx)]))>1
                B(pos(i),:)=-1;
            else
                B(pos(i),:)=unique([Amat(pos(i),AdelIdx) Bmat(pos(i),BdelIdx)]);
            end
            Bcounts(pos(i),:)=AcountsOrig(pos(i),:).*AdelIdx+BcountsOrig(pos(i),:).*BdelIdx;
            BpopAF(pos(i),:)=ApopAForig(pos(i),:).*AdelIdx+BpopAForig(pos(i),:).*BdelIdx;
        elseif(tIdx(2)==3)
            if length(unique([Amat(pos(i),AinsIdx) Bmat(pos(i),BinsIdx)]))>1
                B(pos(i),:)=-1;
            else
                B(pos(i),:)=unique([Amat(pos(i),AinsIdx) Bmat(pos(i),BinsIdx)]);
            end
            Bcounts(pos(i),:)=AcountsOrig(pos(i),:).*AinsIdx+BcountsOrig(pos(i),:).*BinsIdx;
            BpopAF(pos(i),:)=ApopAForig(pos(i),:).*AinsIdx+BpopAForig(pos(i),:).*BinsIdx;
        elseif(tIdx(2)==1)
            if length(unique([Amat(pos(i),ArefIdx) Bmat(pos(i),BrefIdx)]))~=1
                B(pos(i),:)=-1;
            else
                B(pos(i),:)=unique([Amat(pos(i),ArefIdx) Bmat(pos(i),BrefIdx)]);
            end
            Bcounts(pos(i),:)=AcountsOrig(pos(i),:).*ArefIdx+BcountsOrig(pos(i),:).*BrefIdx;
            BpopAF(pos(i),:)=ApopAForig(pos(i),:).*ArefIdx+BpopAForig(pos(i),:).*BrefIdx;
        else
            if length(unique([Amat(pos(i),AaltIdx) Bmat(pos(i),BaltIdx)]))>1
                B(pos(i),:)=-1;
            else
                B(pos(i),:)=unique([Amat(pos(i),AaltIdx) Bmat(pos(i),BaltIdx)]);
            end
            Bcounts(pos(i),:)=AcountsOrig(pos(i),:).*AaltIdx+BcountsOrig(pos(i),:).*BaltIdx;
            BpopAF(pos(i),:)=ApopAForig(pos(i),:).*AaltIdx+BpopAForig(pos(i),:).*BaltIdx;
        end
        if max(tIdx(1:2))==4
            if length(unique(RefOrig(pos(i),AdelIdx | BdelIdx)))>1
                Ref(pos(i),:)=-1;
            else
                Ref(pos(i),:)=unique(RefOrig(pos(i),AdelIdx | BdelIdx));
                if tIdx(1)==1 && tIdx(2)==3
                    B(pos(i),:)=unique(RefOrig(pos(i),AdelIdx | BdelIdx));
                elseif tIdx(2)==1 && tIdx(1)==3
                    A(pos(i),:)=unique(RefOrig(pos(i),AdelIdx | BdelIdx));
                end
            end
        end
    else
        A(pos(i))=alleles(idx(1));
        B(pos(i))=alleles(idx(2));
        for j=1:size(Amat,2)
            if A(pos(i))==Amat(pos(i),j)
                Acounts(pos(i),j)=AcountsOrig(pos(i),j);
                ApopAF(pos(i),j)=ApopAForig(pos(i),j);
            elseif A(pos(i))==Bmat(pos(i),j)
                Acounts(pos(i),j)=BcountsOrig(pos(i),j);
                ApopAF(pos(i),j)=BpopAForig(pos(i),j);
            elseif A(pos(i))==RefOrig(pos(i),j)
                Acounts(pos(i),j)=RDmat(pos(i),j);
                ApopAF(pos(i),j)=max(ApopAF(:,j));
            else
                Acounts(pos(i),j)=0;
                if A(pos(i))>4 | B(pos(i))>4 | RefOrig(pos(i),j)
                    ApopAF(pos(i),j)=inputParam.pvFreqIndel;
                else
                    ApopAF(pos(i),j)=inputParam.pvFreq;
                end
            end
            if A(pos(i))==Amat(pos(i),j) && B(pos(i))==Bmat(pos(i),j)
                Ref(pos(i),j)=RefOrig(pos(i),j);
            end
        end
        for j=1:size(Bmat,2)
            if B(pos(i))==Amat(pos(i),j)
                Bcounts(pos(i),j)=AcountsOrig(pos(i),j);
                BpopAF(pos(i),j)=ApopAForig(pos(i),j);
            elseif B(pos(i))==Bmat(pos(i),j)
                Bcounts(pos(i),j)=BcountsOrig(pos(i),j);
                BpopAF(pos(i),j)=BpopAForig(pos(i),j);
            elseif B(pos(i))==RefOrig(pos(i),j)
                Bcounts(pos(i),j)=RDmat(pos(i),j);
                BpopAF(pos(i),j)=max(BpopAForig(:,j));
            else
                Bcounts(pos(i),j)=0;
                if A(pos(i))>4 | B(pos(i))>4 | RefOrig(pos(i),j)
                    BpopAF(pos(i),j)=inputParam.pvFreqIndel;
                else
                    BpopAF(pos(i),j)=inputParam.pvFreq;
                end
            end
        end
    end
    i;
end

Ref=max(Ref,[],2);
ApopAF=max(ApopAF,[],2);
BpopAF=max(BpopAF,[],2);
cosmic=max(cosmic,[],2);

%%% find priors
indelPos=A>4 | B>4 | Ref >4 | A<0 | B<0 | Ref <0;
priorNonDip=max(mapQC,[],2)-inputParam.minLik;
priorSomatic(indelPos,:)=((cosmic(indelPos)+1).*inputParam.priorSomaticIndel).*(1-priorNonDip(indelPos));
priorSomatic(~indelPos,:)=((cosmic(~indelPos)+1).*inputParam.priorSomaticSNV).*(1-priorNonDip(~indelPos));
priorHet=2*ApopAF.*BpopAF.*(1-priorSomatic-priorNonDip);
priorHom=(ApopAF.^2).*(1-priorSomatic-priorNonDip);
priorOther=1-priorHet-priorHom-priorSomatic-priorNonDip;

prior=array2table([priorSomatic priorHet priorHom priorOther priorNonDip],'VariableNames',{'Somatic','Het','Hom','Other','nonDiploid'});

%%% calculate expected heterozygous allele frequencies
cnCorr=cnaF.*M./N+(1-cnaF)*0.5;
cnCorr(N==0)=0.5;

for i=1:size(Amat,2)
    pos=find(RDmat(:,i)==0);
    idx=getPosInRegions(posList(pos,:),exonRD{i}(:,1:3));
    RDmat(pos(~isnan(idx)),i)=exonRD{i}(idx(~isnan(idx)),4);
end

Acounts(Acounts+Bcounts==0 & (A==Ref)*ones(1,size(Amat,2)))=RDmat(Acounts+Bcounts==0 & (A==Ref)*ones(1,size(Amat,2)));
%%% calculate likliehoods of germline genotypes
for i=1:size(Acounts,2)
    pDataHom(:,i)=bbinopdf_ln(Acounts(:,i),RDmat(:,i),Wmat(:,i).*(1-10.^(-BmeanBQ(:,i)./10)),Wmat(:,i).*10.^(-BmeanBQ(:,i)./10));
    pDataHet(:,i)=bbinopdf_ln(Bcounts(:,i),RDmat(:,i),cnCorr(:,i).*Wmat(:,i),(1-cnCorr(:,i)).*Wmat(:,i));
    pDataOther(:,i)=bbinopdf_ln(max(RDmat(:,i)-Acounts(:,i)-Bcounts(:,i),0),RDmat(:,i),Wmat(:,i).*(1-10.^(-AmeanBQ(:,i)./10)),Wmat(:,i).*10.^(-AmeanBQ(:,i)./10));
end

%%% calculate expected somatic allele frequencies
for j=1:size(Acounts,2)
    for i=1:length(f(j,:))
        matchIdx=round(100.*cnaF(:,j))==round(100.*f(j,i));
        expAF(matchIdx,i,j)=f(j,i)*(N(matchIdx,j)-M(matchIdx,j))./(f(j,i)*N(matchIdx,j)+(1-f(j,i))*2);
        if sum(~matchIdx)>0
            expAF(~matchIdx,i,j)=f(j,i)./(cnaF(~matchIdx,j).*N(~matchIdx,j)+(1-cnaF(~matchIdx,j))*2);
            expAF(N(:,j)==0 & ~matchIdx,i,j)=min([(1-cnaF(:,j)); f(j,i)])./2;
        end
    end
end

%%% find likelihood of somatic mutations
tumorIdx=setdiff([1:size(Acounts,2)],inputParam.NormalSample);
bIdx=ApopAF>=BpopAF;
aIdx=ApopAF<BpopAF;
meanAF=mean(Bcounts(:,tumorIdx)./RDmat(:,tumorIdx),2);
for j=1:size(Acounts,2)
    if j==inputParam.NormalSample
        cloneLik(:,:,j)=pDataHom(:,j)*ones(1,size(Acounts,2));
    else
        for i=1:length(f(j,:))
            alpha(:,i)=max(expAF(:,i,j),inputParam.minLik)*W(j);
            beta(:,i)=(1-max(expAF(:,i,j),inputParam.minLik))*W(j);
            cloneLik(bIdx,i,j)=bbinopdf_ln(Bcounts(bIdx,j),RDmat(bIdx,j),alpha(bIdx,i),beta(bIdx,i));
            cloneLik(aIdx,i,j)=bbinopdf_ln(Acounts(aIdx,j),RDmat(aIdx,j),alpha(aIdx,i),beta(aIdx,i));
        end
    end
    if inputParam.NormalSample>0
        pDataNonDip(:,j)=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),meanAF.*W(j),(1-meanAF).*W(j));
    else
        pDataNonDip(:,j)=zeros(size(Bcounts(:,j)));
    end
end
[~,cloneId]=max(prod(cloneLik,3),[],2);
for i=1:size(f,2)
    pDataSomatic(cloneId==i,:)=squeeze(cloneLik(cloneId==i,i,:));
end
pDataNonDip(isnan(pDataNonDip))=0;



%condP=array2table([pDataSomatic pDataHet pDataHom pDataOther],'VariableNames',{'Somatic','Het','Hom','Other'});


%%% find expected somatic AF for most likely clone
% for i=1:length(f)
%     expectedAF(cloneId==i,:)=expAF(cloneId==i,i);
% end

%%% calculate posteriors


pData=prod(pDataHom,2).*priorHom+prod(pDataHet,2).*priorHet+prod(pDataSomatic,2).*priorSomatic+prod(pDataOther,2).*priorOther+prod(pDataNonDip,2).*priorNonDip;
pSomatic=prod(pDataSomatic,2).*priorSomatic./pData;
pGermline=prod(pDataHet,2).*priorHet./pData;
pHom=prod(pDataHom,2).*priorHom./pData;
pOther=prod(pDataOther,2).*priorOther./pData;
pNonDip=prod(pDataNonDip,2).*priorNonDip./pData;

postComb=array2table([posList pSomatic pGermline pHom pOther pNonDip],'VariableNames',{'Chr', 'Pos', 'Somatic','Het','Hom','Other','NonDip'});
countsAll=array2table([posList Ref A B],'VariableNames',{'Chr','Pos','Ref','A','B'});
for i=1:size(Acounts,2)
    pDataComb{i}=array2table([posList pDataSomatic(:,i) pDataHet(:,i) pDataHom(:,i) pDataOther(:,i) pDataNonDip(:,i)],'VariableNames',{'Chr', 'Pos', 'Somatic','Het','Hom','Other','NonDip'});
    countsAll.Acounts(:,i)=Acounts(:,i);
    countsAll.Bcounts(:,i)=Bcounts(:,i);
end

pDataSum=sum(pData);
    

            