function [postComb, pDataSum, pDataComb,cloneId,prior,alleleId]=jointSNV(Tcell, fIn, W, inputParam)
%jointSNV - find joint probability that variants are somatic, germline
%heterozygous or homozygous across samples
% Syntax:  [postComb, pDataSum, pDataComb,cloneId,prior,alleleId]=jointSNV(Tcell, exonRD, fIn, W, inputParam)
%
% Inputs:
%    Tcell - cell array of tables with length equal to the number of bams,
%       each table must have the following columns:  {'Chr','Pos', 'ReadDepthPass',
%       'RefComb','AComb','AcountsComb','AmeanBQ','BComb','BCountsComb','BmeanBQ',
%       'ApopAFcomb','BpopAFcomb','CosmicCount','PosMapQC'};
%    fIn - matrix of sample fractions with one row for each tumor sample in
%          Tcell and one column for each clonal variant group
%    W - parameter controling spread of read count distributions
%    inputParam - data structure with fields from paramFile
%
% Outputs:
%    postComb - table of posterior probabilites that variants belong to
%               each class
%    pDataSum - sum of probability of data
%    pDataComb - table of probability of data given each variant class
%                 model
%    cloneId - index of column of fIn that best fits data given that the
%              variant is somatic
%    prior - table of prior probabilities that variant belongs to each
%            class
%    alleleId - indicates whether the variant is on the minor (1) or 
%               major (2) allele given the variant is somatic
%
% Other m-files required: none
% Other requirements: none
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain, getCounts, preprocess?

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 7-May-2018

%------------- BEGIN CODE --------------
%%% fill in values of F for normal sample
f=zeros(length(Tcell),inputParam.numClones);
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
f(tIdx,1:end)=reshape(fIn,[],inputParam.numClones);

%%%% put data from Tcell into matrices
A=NaN(size(Tcell{1},1),1);
B=NaN(size(Tcell{1},1),1);
Acounts=NaN(size(Tcell{1},1),length(Tcell));
Bcounts=NaN(size(Tcell{1},1),length(Tcell));
ApopAF=NaN(size(Tcell{1},1),1);
BpopAF=NaN(size(Tcell{1},1),1);
RDmat=NaN(size(Tcell{1},1),length(Tcell));
Ref=NaN(size(Tcell{1},1),1);
cosmic=NaN(size(Tcell{1},1),length(Tcell));
cnaF=NaN(size(Tcell{1},1),length(Tcell));
N=NaN(size(Tcell{1},1),1);
M=NaN(size(Tcell{1},1),1);
Wmat=NaN(size(Tcell{1},1),length(Tcell));
BmeanBQ=NaN(size(Tcell{1},1),length(Tcell));
AmeanBQ=NaN(size(Tcell{1},1),length(Tcell));
mapQC=NaN(size(Tcell{1},1),length(Tcell));
for i=1:size(Tcell,2)
    T=Tcell{i};
    %[lia,locb]=ismember([T.Chr T.Pos],posList,'rows');
    locb=1:height(T);
    A(locb,1)=T.Acomb;
    B(locb,1)=T.Bcomb;
    Acounts(locb,i)=T.AcountsComb;
    Bcounts(locb,i)=T.BcountsComb;
    ApopAF(locb,1)=T.ApopAFcomb;
    BpopAF(locb,1)=T.BpopAFcomb;
    RDmat(locb,i)=T.ReadDepthPass;
    Ref(locb,1)=T.RefComb;
    cosmic(locb,i)=T.CosmicCount;
    cnaF(locb,i)=T.cnaF;
    N(locb,i)=T.NumCopies;
    M(locb,i)=T.MinAlCopies;
    Wmat(locb,i)=T.W;
    BmeanBQ(locb,i)=T.BmeanBQ;
    AmeanBQ(locb,i)=T.AmeanBQ;
    mapQC(locb,i)=T.PosMapQC;
end
mapQC(isnan(mapQC))=10.^(inputParam.minPosQual./-10);
M=min(M,N-M);
cosmic=max(cosmic,[],2);
for i=1:size(Wmat,2)
    Wmat(Wmat(:,i)==0,i)=mode(Wmat(Wmat(:,i)~=0,i));
    BmeanBQ(BmeanBQ(:,i)==0,i)=inputParam.defaultBQ;
    AmeanBQ(AmeanBQ(:,i)==0,i)=inputParam.defaultBQ;
end

%%% find priors
idx=Acounts+Bcounts>RDmat;
RDmat(idx)=Acounts(idx)+Bcounts(idx);
germAF=NaN(size(Acounts,1),1);
if inputParam.NormalSample>0
    bIdx=Acounts(:,inputParam.NormalSample)>=Bcounts(:,inputParam.NormalSample);
    aIdx=Acounts(:,inputParam.NormalSample)<Bcounts(:,inputParam.NormalSample);
    germAF(bIdx,:)=Bcounts(bIdx,inputParam.NormalSample)./(RDmat(bIdx,inputParam.NormalSample));
    germAF(aIdx,:)=Acounts(aIdx,inputParam.NormalSample)./(RDmat(aIdx,inputParam.NormalSample));
    germAF(RDmat(:,inputParam.NormalSample)==0,:)=0;
    tAF=mean(Bcounts(:,tIdx),2)./mean(RDmat(:,tIdx),2);
    germAF=min(germAF,tAF);
else
    bIdx=ApopAF>=BpopAF;
    aIdx=ApopAF<BpopAF;
end
indelPos=A>4 | B>4 | Ref >4 | A<0 | B<0 | Ref <0;
priorNonDip=min(max(max(mapQC,[],2),inputParam.minLik),0.5);
priorHom=(ApopAF.^2).*(1-priorNonDip);
priorSomatic(indelPos & aIdx,:)=((BpopAF(indelPos & aIdx).^2).*(cosmic(indelPos & aIdx)+1).*inputParam.priorSomaticIndel).*(1-priorNonDip(indelPos & aIdx));
priorSomatic(indelPos & bIdx,:)=((ApopAF(indelPos & bIdx).^2).*(cosmic(indelPos & bIdx)+1).*inputParam.priorSomaticIndel).*(1-priorNonDip(indelPos & bIdx));
priorSomatic(~indelPos & aIdx,:)=((BpopAF(~indelPos & aIdx).^2).*(cosmic(~indelPos & aIdx)+1).*inputParam.priorSomaticSNV).*(1-priorNonDip(~indelPos & aIdx));
priorSomatic(~indelPos & bIdx,:)=((ApopAF(~indelPos & bIdx).^2).*(cosmic(~indelPos & bIdx)+1).*inputParam.priorSomaticSNV).*(1-priorNonDip(~indelPos & bIdx));
priorHet=2*ApopAF.*BpopAF.*(1-priorSomatic-priorNonDip);
priorOther=1-priorHet-priorHom-priorSomatic-priorNonDip;
prior=array2table([priorSomatic priorHet priorHom priorOther priorNonDip],'VariableNames',{'Somatic','Het','Hom','Other','nonDiploid'});

%%% calculate expected heterozygous allele frequencies
expAFhet=(cnaF.*M+(1-cnaF))./(cnaF.*N+(1-cnaF)*2);
expAFhet(N==0 & cnaF<1)=0.5;
expAFhet(N==0 & cnaF==1)=0;
expAFhet=max(expAFhet,10.^(inputParam.defaultBQ./-10));

%%% calculate likliehoods of germline genotypes
pDataHom=NaN(size(Acounts));
pDataHet=NaN(size(Acounts));
pDataOther=NaN(size(Acounts));
for i=1:size(Acounts,2)
    pDataHom(:,i)=bbinopdf_ln(Acounts(:,i),RDmat(:,i),Wmat(:,i).*(1-10.^(-AmeanBQ(:,i)./10)),Wmat(:,i).*10.^(-AmeanBQ(:,i)./10));
    pDataHet(:,i)=bbinopdf_ln(Bcounts(:,i),RDmat(:,i),expAFhet(:,i).*Wmat(:,i),(1-expAFhet(:,i)).*Wmat(:,i))+(binopdf(round(expAFhet(:,i).*Wmat(:,i)),round(Wmat(:,i)),1-expAFhet(:,i))./binopdf(round(expAFhet(:,i).*Wmat(:,i)),round(Wmat(:,i)),expAFhet(:,i))).*bbinopdf_ln(Acounts(:,i),RDmat(:,i),expAFhet(:,i).*Wmat(:,i),(1-expAFhet(:,i)).*Wmat(:,i));
    pDataOther(:,i)=bbinopdf_ln(max(RDmat(:,i)-Acounts(:,i)-Bcounts(:,i),0),RDmat(:,i),Wmat(:,i).*(1-10.^(-AmeanBQ(:,i)./10)),Wmat(:,i).*10.^(-AmeanBQ(:,i)./10));
end

%%% calculate expected somatic allele frequencies
expAF=NaN(size(Acounts,1),size(f,2),size(Acounts,2),2);
for j=1:size(Acounts,2)
    for i=1:length(f(j,:))
        matchIdx=round(100.*cnaF(:,j))==round(100.*f(j,i));
        expAF(matchIdx,i,j,1)=f(j,i)*(M(matchIdx,j))./(f(j,i)*N(matchIdx,j)+(1-f(j,i))*2);
        expAF(matchIdx,i,j,2)=f(j,i)*(N(matchIdx,j)-M(matchIdx,j))./(f(j,i)*N(matchIdx,j)+(1-f(j,i))*2);
        if sum(~matchIdx)>0
            expAF(~matchIdx,i,j,1)=f(j,i)./(cnaF(~matchIdx,j).*N(~matchIdx,j)+(1-cnaF(~matchIdx,j))*2);
       end
    end
end
expAF(expAF>1)=1;
expAF(expAF<0)=0;

%%% find likelihood of somatic mutations
cloneLik=NaN(size(expAF));
alpha=NaN(size(Acounts,1),size(f,2));
beta=NaN(size(Acounts,1),size(f,2));
shiftAF=NaN(size(Acounts,1),size(Acounts,2),2);
pDataNonDip=NaN([size(Acounts) 2]);
for j=1:size(Acounts,2)
    if j==inputParam.NormalSample
        cloneLik(bIdx,:,j,1)=bbinopdf_ln(Acounts(bIdx,j),RDmat(bIdx,j),Wmat(bIdx,j).*(1-10.^(-AmeanBQ(bIdx,j)./10)),Wmat(bIdx,j).*10.^(-AmeanBQ(bIdx,j)./10))*ones(1,length(f(j,:)));
        cloneLik(aIdx,:,j,1)=bbinopdf_ln(Bcounts(aIdx,j),RDmat(aIdx,j),Wmat(aIdx,j).*(1-10.^(-BmeanBQ(aIdx,j)./10)),Wmat(aIdx,j).*10.^(-BmeanBQ(aIdx,j)./10))*ones(1,length(f(j,:)));
        cloneLik(bIdx,:,j,2)=bbinopdf_ln(Acounts(bIdx,j),RDmat(bIdx,j),Wmat(bIdx,j).*(1-10.^(-AmeanBQ(bIdx,j)./10)),Wmat(bIdx,j).*10.^(-AmeanBQ(bIdx,j)./10))*ones(1,length(f(j,:)));
        cloneLik(aIdx,:,j,2)=bbinopdf_ln(Bcounts(aIdx,j),RDmat(aIdx,j),Wmat(aIdx,j).*(1-10.^(-BmeanBQ(aIdx,j)./10)),Wmat(aIdx,j).*10.^(-BmeanBQ(aIdx,j)./10))*ones(1,length(f(j,:)));
    else
        for i=1:length(f(j,:))
            for k=1:2
                alpha(:,i)=max(expAF(:,i,j,k),inputParam.minLik)*W(j);
                beta(:,i)=(1-max(expAF(:,i,j,k),inputParam.minLik))*W(j);
                cloneLik(bIdx,i,j,k)=bbinopdf_ln(Bcounts(bIdx,j),RDmat(bIdx,j),alpha(bIdx,i),beta(bIdx,i));
                cloneLik(aIdx,i,j,k)=bbinopdf_ln(Acounts(aIdx,j),RDmat(aIdx,j),alpha(aIdx,i),beta(aIdx,i));
            end
        end
    end
    if inputParam.NormalSample>0
        shiftAF(:,j,1)=(cnaF(:,inputParam.NormalSample).*N(:,inputParam.NormalSample)./M(:,inputParam.NormalSample)+2*(1-cnaF(:,inputParam.NormalSample))).*cnaF(:,j).*(M(:,j)./N(:,j)).*germAF+(1-cnaF(:,j)).*germAF;
        shiftAF(M(:,inputParam.NormalSample)==0,j,1)=0;
        shiftAF(:,j,2)=(cnaF(:,inputParam.NormalSample).*N(:,inputParam.NormalSample)./(N(:,inputParam.NormalSample)-M(:,inputParam.NormalSample))+2*(1-cnaF(:,inputParam.NormalSample))).*cnaF(:,j).*((N(:,j)-M(:,j))./N(:,j)).*germAF+(1-cnaF(:,j)).*germAF;
        shiftAF(~isfinite(shiftAF(:,j,1)),j,1)=germAF(~isfinite(shiftAF(:,j,1)));
        shiftAF(~isfinite(shiftAF(:,j,2)),j,2)=germAF(~isfinite(shiftAF(:,j,2)));
        pDataNonDip(:,j,1)=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),min(shiftAF(:,j,1),1-shiftAF(:,j,1)).*W(j),max(shiftAF(:,j,1),1-shiftAF(:,j,1)).*W(j));
        pDataNonDip(:,j,2)=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),min(shiftAF(:,j,2),1-shiftAF(:,j,2)).*W(j),max(shiftAF(:,j,2),1-shiftAF(:,j,2)).*W(j));
    else
        pDataNonDip(:,j,1)=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),germAF.*W(j),(1-germAF).*W(j));
        pDataNonDip(:,j,2)=pDataNonDip(:,j,1);
    end
end

%%%find most likely clone and allele
cloneLik(isnan(expAF))=NaN;
[~,nonDipIdx]=max(prod(pDataNonDip,2),[],3);
pDataNonDipMax=NaN(size(Acounts));
for i=1:size(pDataNonDip,3)
    pDataNonDipMax(nonDipIdx==i,:)=pDataNonDip(nonDipIdx==i,:,i);
end
[m,mIdx]=max(prod(cloneLik,3),[],4);
[~,cloneId]=max(m,[],2);
alleleId=NaN(size(cloneId));
for i=1:size(cloneLik,2)
    alleleId(cloneId==i,:)=mIdx(cloneId==i,i);
end

%%%%find maximum somatic likelihood
pDataSomatic=NaN(size(Acounts));
expAFsom=NaN(size(Acounts));
for i=1:size(f,2)
    for k=1:2
        pDataSomatic(cloneId==i & alleleId==k,:)=squeeze(cloneLik(cloneId==i & alleleId==k,i,:,k));
        expAFsom(cloneId==i & alleleId==k,:)=squeeze(expAF(cloneId==i & alleleId==k,i,:,k));
    end
end
pDataSomatic(Bcounts<1)=min(pDataSomatic(Bcounts<1),pDataHom(Bcounts<1));
pDataNonDipMax(isnan(pDataNonDipMax))=0;

%%%% find posterior probabilites
pData=prod(pDataHom,2).*priorHom+prod(pDataHet,2).*priorHet+prod(pDataSomatic,2).*priorSomatic+prod(pDataOther,2).*priorOther+prod(pDataNonDipMax,2).*priorNonDip;
pSomatic=prod(pDataSomatic,2).*priorSomatic./pData;
pGermline=prod(pDataHet,2).*priorHet./pData;
pHom=prod(pDataHom,2).*priorHom./pData;
pOther=prod(pDataOther,2).*priorOther./pData;
pNonDip=prod(pDataNonDipMax,2).*priorNonDip./pData;
postComb=array2table([Tcell{1}.Chr Tcell{1}.Pos pSomatic pGermline pHom pOther pNonDip],'VariableNames',{'Chr', 'Pos', 'Somatic','Het','Hom','Other','NonDip'});
pDataComb=cell(size(Tcell));
for i=1:size(Acounts,2)
    pDataComb{i}=array2table([Tcell{1}.Chr Tcell{1}.Pos pDataSomatic(:,i) pDataHet(:,i) pDataHom(:,i) pDataOther(:,i) pDataNonDipMax(:,i)],'VariableNames',{'Chr', 'Pos', 'Somatic','Het','Hom','Other','NonDip'});
end
pDataSum=sum(pData);
    

            
