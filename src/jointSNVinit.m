function [postComb, pDataSum, somaticFlag, pDataComb,cloneId,prior]=jointSNVinit(Tcell, exonRD, f, W, inputParam)

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
B(min(Bmat,[],2)-max(Bmat,[],2)==0,:)=Amat(min(Bmat,[],2)-max(Bmat,[],2)==0,1);
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
    A(pos(i))=alleles(idx==1);
    B(pos(i))=alleles(idx==2);
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
    i;
end

Ref=max(Ref,[],2);
ApopAF=max(ApopAF,[],2);
BpopAF=max(BpopAF,[],2);
cosmic=max(cosmic,[],2);

%%% find priors
indelPos=A>4 | B>4 | Ref >4;
priorNonDip=max(mapQC,[],2);
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
    for i=1:length(f(j,:))
        alpha(:,i)=expAF(:,i,j)*W(j);
        beta(:,i)=(1-expAF(:,i,j))*W(j);
        cloneLik(bIdx,i)=bbinopdf_ln(Bcounts(bIdx,j),RDmat(bIdx,j),alpha(bIdx,i),beta(bIdx,i));
        cloneLik(aIdx,i)=bbinopdf_ln(Acounts(aIdx,j),RDmat(aIdx,j),alpha(aIdx,i),beta(aIdx,i));
    end
    [pDataSomatic(:,j),cloneId(:,j)]=max(cloneLik,[],2);
    if inputParam.NormalSample>0
        pDataNonDip(:,j)=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),meanAF.*W(j),(1-meanAF).*W(j));
    else
        pDataNonDip(:,j)=zeros(size(alpha2));
    end
end
pDataNonDip(isnan(pDataNonDip))=0;



%condP=array2table([pDataSomatic pDataHet pDataHom pDataOther],'VariableNames',{'Somatic','Het','Hom','Other'});


%%% find expected somatic AF for most likely clone
% for i=1:length(f)
%     expectedAF(cloneId==i,:)=expAF(cloneId==i,i);
% end

%%% calculate posteriors

[m,idx]=max(pDataSomatic(:,tumorIdx),[],2);
for i=1:size(pDataSomatic,2)
    if inputParam.NormalSample==i
        pDataSomaticComb(:,i)=pDataHom(:,i);
    else
        pDataSomaticComb(i==tumorIdx(idx),i)=pDataSomatic(i==tumorIdx(idx),i);
        pDataSomaticComb(i~=tumorIdx(idx),i)=max(pDataSomatic(i~=tumorIdx(idx),i),pDataHom(i~=tumorIdx(idx),i));
    end
end
somaticFlag=zeros(size(pDataSomaticComb));
somaticFlag(pDataSomatic==pDataSomaticComb)=1;


['priorHom: ' num2str(size(priorHom))];
['pDataHom: ' num2str(size(pDataHom))];
['priorHet: ' num2str(size(priorHet))];
['pDataHet: ' num2str(size(pDataHet))];
['pDataSomaticComb: ' num2str(size(pDataSomaticComb))];
['priorSomatic: ' num2str(size(priorSomatic))];
['pDataOther: ' num2str(size(pDataOther))];
['priorOther: ' num2str(size(priorOther))];
pData=prod(pDataHom,2).*priorHom+prod(pDataHet,2).*priorHet+prod(pDataSomaticComb,2).*priorSomatic+prod(pDataOther,2).*priorOther+prod(pDataNonDip,2).*priorNonDip;
pSomatic=prod(pDataSomaticComb,2).*priorSomatic./pData;
pGermline=prod(pDataHet,2).*priorHet./pData;
pHom=prod(pDataHom,2).*priorHom./pData;
pOther=prod(pDataOther,2).*priorOther./pData;
pNonDip=prod(pDataNonDip,2).*priorNonDip./pData;

postComb=array2table([posList pSomatic pGermline pHom pOther pNonDip],'VariableNames',{'Chr', 'Pos', 'Somatic','Het','Hom','Other','NonDip'});
for i=1:size(Acounts,2)
    pDataComb{i}=array2table([posList pDataSomatic(:,i) pDataHet(:,i) pDataHom(:,i) pDataOther(:,i) pDataNonDip(:,i)],'VariableNames',{'Chr', 'Pos', 'Somatic','Het','Hom','Other','NonDip'});
end

pDataSum=sum(pData);
    

            