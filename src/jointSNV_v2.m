function [postComb, pDataSum, pDataComb,cloneId,prior,countsAll]=jointSNV_v2(Tcell, exonRD, fIn, W, inputParam)

f=zeros(length(Tcell),inputParam.numClones);
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
f(tIdx,1:end)=reshape(fIn,[],inputParam.numClones);

posList=[];
for i=1:size(Tcell,2)
    T=Tcell{i};
    posList=[posList; T.Chr T.Pos];
end

posList=unique(posList,'rows');

for i=1:size(Tcell,2)
    T=Tcell{i};
    [lia,locb]=ismember([T.Chr T.Pos],posList,'rows');
    A(locb,1)=T.Acomb;
    B(locb,1)=T.Bcomb;
    Acounts(locb,i)=T.AcountsComb;
    Bcounts(locb,i)=T.BcountsComb;
    ApopAF(locb,1)=T.ApopAFcomb;
    BpopAF(locb,1)=T.BpopAFcomb;
    RDmat(locb,i)=T.ReadDepthPass;
    Ref(locb,1)=T.RefComb;
    cosmic(locb,1)=T.CosmicCount;
    cnaF(locb,i)=T.cnaF;
    N(locb,i)=T.NumCopies;
    M(locb,i)=T.MinAlCopies;
    Wmat(locb,i)=T.W;
    BmeanBQ(locb,i)=T.BmeanBQ;
    AmeanBQ(locb,i)=T.AmeanBQ;
    mapQC(locb,i)=T.PosMapQC;
end

for i=1:size(Wmat,2)
    Wmat(Wmat(:,i)==0,i)=mode(Wmat(Wmat(:,i)~=0,i));
    BmeanBQ(BmeanBQ(:,i)==0,i)=inputParam.defaultBQ;
    AmeanBQ(AmeanBQ(:,i)==0,i)=inputParam.defaultBQ;
end

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
cnCorr=max(cnCorr,10.^(inputParam.defaultBQ./-10));

for i=1:size(Acounts,2)
    pos=find(RDmat(:,i)==0);
    idx=getPosInRegions(posList(pos,:),exonRD{i}(:,1:3));
    RDmat(pos(~isnan(idx)),i)=exonRD{i}(idx(~isnan(idx)),4);
end

Acounts(Acounts+Bcounts==0 & (A==Ref)*ones(1,size(Acounts,2)))=RDmat(Acounts+Bcounts==0 & (A==Ref)*ones(1,size(Acounts,2)));
%%% calculate likliehoods of germline genotypes
for i=1:size(Acounts,2)
    pDataHom(:,i)=bbinopdf_ln(Acounts(:,i),RDmat(:,i),Wmat(:,i).*(1-10.^(-AmeanBQ(:,i)./10)),Wmat(:,i).*10.^(-AmeanBQ(:,i)./10));
    pDataHet(:,i)=bbinopdf_ln(Bcounts(:,i),RDmat(:,i),cnCorr(:,i).*Wmat(:,i),(1-cnCorr(:,i)).*Wmat(:,i));
    pDataOther(:,i)=bbinopdf_ln(max(RDmat(:,i)-Acounts(:,i)-Bcounts(:,i),0),RDmat(:,i),Wmat(:,i).*(1-10.^(-AmeanBQ(:,i)./10)),Wmat(:,i).*10.^(-AmeanBQ(:,i)./10));
end

%%% calculate expected somatic allele frequencies
for j=1:size(Acounts,2)
    for i=1:length(f(j,:))
        matchIdx=round(100.*cnaF(:,j))==round(100.*f(j,i));
        subIdx=cnaF(:,j)+f(j,i)>1 & f(j,i)<cnaF(:,j);
        expAF(matchIdx,i,j)=f(j,i)*(N(matchIdx,j)-M(matchIdx,j))./(f(j,i)*N(matchIdx,j)+(1-f(j,i))*2);
        if sum(~matchIdx)>0
            expAF(~matchIdx,i,j)=f(j,i)./(cnaF(~matchIdx,j).*N(~matchIdx,j)+(1-cnaF(~matchIdx,j))*2);
            expAF(N(:,j)==0 & ~matchIdx,i,j)=f(j,i)./2;
            expAF(subIdx,i,j)=f(j,i)*M(subIdx,j)./(cnaF(subIdx,j).*N(subIdx,j)+(1-cnaF(subIdx,j))*2);
            %expAF(N(:,j)==0 & ~matchIdx,i,j)=min((1-cnaF(N(:,j)==0 & ~matchIdx,j)),f(j,i))./2;
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
        cloneLik(:,:,j)=pDataHom(:,j)*ones(1,length(f(j,:)));
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
pDataSomatic(Bcounts<1)=min(pDataSomatic(Bcounts<1),pDataHom(Bcounts<1));
pDataNonDip(isnan(pDataNonDip))=0;
if inputParam.NormalSample>0
    pDataNonDip(Bcounts(:,inputParam.NormalSample)==0,inputParam.NormalSample)=0;
end



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
    

            