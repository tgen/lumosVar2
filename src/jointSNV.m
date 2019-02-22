function [postComb, pDataSum, cloneId,prior,alleleId,pDataSomatic]=jointSNV(posTable, fIn, W, inputParam)
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
f=zeros(inputParam.sampleCount,inputParam.numClones);
tIdx=setdiff(1:inputParam.sampleCount,inputParam.NormalSample);
f(tIdx,1:end)=reshape(fIn,[],inputParam.numClones);


%%%% put data from Tcell into matrices
A=posTable.A_1;
B=posTable.B_1;
Acounts=posTable{:,strncmp(posTable.Properties.VariableNames,'AcountsComb_',12)};
Bcounts=posTable{:,strncmp(posTable.Properties.VariableNames,'BcountsComb_',12)};
ApopAF=max(posTable.ApopAF_1,inputParam.pvFreqIndel);
BpopAF=max(posTable.BpopAF_1,inputParam.pvFreqIndel);
RDmat=posTable{:,strncmp(posTable.Properties.VariableNames,'ReadDepthPass_',14)};
Ref=posTable.Ref_1;
cosmic=posTable.CosmicCount_1;
cnaF=posTable{:,strncmp(posTable.Properties.VariableNames,'cnaF_',5)};
N=posTable.N;
M=posTable.M;
Wmat=posTable{:,strncmp(posTable.Properties.VariableNames,'W_',2)};
BmeanBQ=posTable{:,strncmp(posTable.Properties.VariableNames,'BmeanBQ_',8)};
AmeanBQ=posTable{:,strncmp(posTable.Properties.VariableNames,'AmeanBQ_',8)};
mapQC=posTable.PosMapQC_1;
aIdx=posTable.aIdx;
bIdx=~posTable.aIdx;

mapQC(isnan(mapQC))=10.^(inputParam.minPosQual./-10);
M=min(M,N-M);
cosmic=max(cosmic,[],2);
for i=1:inputParam.sampleCount
    %Wmat(Wmat(:,i)==0,i)=mode(Wmat(Wmat(:,i)~=0,i));
    BmeanBQ(BmeanBQ(:,i)==0,i)=inputParam.defaultBQ;
    AmeanBQ(AmeanBQ(:,i)==0,i)=inputParam.defaultBQ;
    idx=Acounts(:,i)+Bcounts(:,i)>RDmat(:,i);
    RDmat(idx,i)=Acounts(idx,i)+Bcounts(idx,i);
end

%%% find priors
%germAF=NaN(size(Acounts,1),1);
if inputParam.NormalSample>0
    germAF=Bcounts(:,inputParam.NormalSample)./(RDmat(:,inputParam.NormalSample));
    germAF(aIdx,:)=Acounts(aIdx,inputParam.NormalSample)./(RDmat(aIdx,inputParam.NormalSample));
    germAF(RDmat(:,inputParam.NormalSample)==0,:)=0;
    tAF=mean(Bcounts(:,tIdx),2)./mean(RDmat(:,tIdx),2);
    germAF=min(germAF,tAF);
else
    germAF=min(Acounts,0)./0;
end
indelPos=A>4 | B>4 | Ref >4 | A<0 | B<0 | Ref <0;
priorNonDip=min(max(max(mapQC,[],2),inputParam.minLik),0.5);
priorHom=(ApopAF.^2).*(1-priorNonDip);
priorSomatic=((ApopAF.^2).*(cosmic+1).*inputParam.priorSomaticSNV).*(1-priorNonDip);
priorSomatic(~indelPos & aIdx,:)=((BpopAF(~indelPos & aIdx).^2).*(cosmic(~indelPos & aIdx)+1).*inputParam.priorSomaticSNV).*(1-priorNonDip(~indelPos & aIdx));
priorSomatic(indelPos & aIdx,:)=((BpopAF(indelPos & aIdx).^2).*(cosmic(indelPos & aIdx)+1).*inputParam.priorSomaticIndel).*(1-priorNonDip(indelPos & aIdx));
priorSomatic(indelPos & bIdx,:)=((ApopAF(indelPos & bIdx).^2).*(cosmic(indelPos & bIdx)+1).*inputParam.priorSomaticIndel).*(1-priorNonDip(indelPos & bIdx));
priorHet=2*ApopAF.*BpopAF.*(1-priorSomatic-priorNonDip);
priorOther=1-priorHet-priorHom-priorSomatic-priorNonDip;
prior=array2table([priorSomatic priorHet priorHom priorOther priorNonDip],'VariableNames',{'Somatic','Het','Hom','Other','nonDiploid'});

%%% calculate expected heterozygous allele frequencies
expAFhet=(cnaF.*M+(1-cnaF))./(cnaF.*N+(1-cnaF)*2);
expAFhet(N==0 & cnaF<1)=0.5;
expAFhet(N==0 & cnaF==1)=0;
expAFhet=max(expAFhet,10.^(inputParam.defaultBQ./-10));

%%% calculate likliehoods of germline genotypes
pDataHom=min(Acounts,0);
pDataHet=min(Acounts,0);
pDataOther=min(Acounts,0);
for i=1:inputParam.sampleCount
    pDataHom(:,i)=bbinopdf_ln(Acounts(:,i),RDmat(:,i),Wmat(:,i).*(1-10.^(-AmeanBQ(:,i)./10)),Wmat(:,i).*10.^(-AmeanBQ(:,i)./10));
    pDataHet(:,i)=bbinopdf_ln(Bcounts(:,i),RDmat(:,i),expAFhet(:,i).*Wmat(:,i),(1-expAFhet(:,i)).*Wmat(:,i))+(bbinopdf_ln(round(expAFhet(:,i).*Wmat(:,i)),round(Wmat(:,i)),(1-expAFhet(:,i)).*Wmat(:,i),expAFhet(:,i).*Wmat(:,i))./bbinopdf_ln(round(expAFhet(:,i).*Wmat(:,i)),round(Wmat(:,i)),expAFhet(:,i).*Wmat(:,i),(1-expAFhet(:,i)).*Wmat(:,i))).*bbinopdf_ln(Acounts(:,i),RDmat(:,i),expAFhet(:,i).*Wmat(:,i),(1-expAFhet(:,i)).*Wmat(:,i));
    pDataOther(:,i)=bbinopdf_ln(max(RDmat(:,i)-Acounts(:,i)-Bcounts(:,i),0),RDmat(:,i),Wmat(:,i).*(1-10.^(-AmeanBQ(:,i)./10)),Wmat(:,i).*10.^(-AmeanBQ(:,i)./10));
end

%%% calculate expected somatic allele frequencies
expAF=posTable(:,cellstr({"Chr_1","Pos_1"}));
for j=1:inputParam.sampleCount
    for i=1:length(f(j,:))
        matchIdx=round(100.*cnaF(:,j))==round(100.*f(j,i));
        expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A1'))=f(j,i)*M./(f(j,i)*N+(1-f(j,i))*2);
        expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A2'))=f(j,i)*(N-M)./(f(j,i)*N+(1-f(j,i))*2);
        %if sum(~matchIdx)>0
            expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A1'))(~matchIdx)=f(j,i)./(cnaF(~matchIdx,j).*N(~matchIdx)+(1-cnaF(~matchIdx,j))*2);
        %end
       expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A1'))=min(max(expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A1')),0),1);
       expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A2'))=min(max(expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A2')),0),1);
    end
end

%%% find likelihood of somatic mutations

for i=1:size(f,2) 
    for j=1:inputParam.sampleCount
        if j==inputParam.NormalSample
            currCloneLik1=bbinopdf_ln(Acounts(:,j),RDmat(:,j),Wmat(:,j).*(1-10.^(-AmeanBQ(:,j)./10)),Wmat(:,j).*10.^(-AmeanBQ(:,j)./10));
            currCloneLik1(aIdx)=bbinopdf_ln(Bcounts(aIdx,j),RDmat(aIdx,j),Wmat(aIdx,j).*(1-10.^(-BmeanBQ(aIdx,j)./10)),Wmat(aIdx,j).*10.^(-BmeanBQ(aIdx,j)./10));
            currCloneLik2=bbinopdf_ln(Acounts(:,j),RDmat(:,j),Wmat(:,j).*(1-10.^(-AmeanBQ(:,j)./10)),Wmat(:,j).*10.^(-AmeanBQ(:,j)./10));
            currCloneLik2(aIdx)=bbinopdf_ln(Bcounts(aIdx,j),RDmat(aIdx,j),Wmat(aIdx,j).*(1-10.^(-BmeanBQ(aIdx,j)./10)),Wmat(aIdx,j).*10.^(-BmeanBQ(aIdx,j)./10));
        else
            alpha=max(expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A1')),inputParam.minLik)*W(j);
            beta=(1-max(expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A1')),inputParam.minLik))*W(j);
            currCloneLik1=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),alpha,beta);
            currCloneLik1(aIdx)=bbinopdf_ln(Acounts(aIdx,j),RDmat(aIdx,j),alpha(aIdx),beta(aIdx));
            alpha=max(expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A2')),inputParam.minLik)*W(j);
            beta=(1-max(expAF.(strcat('S',num2str(j),'_C',num2str(i),'_A2')),inputParam.minLik))*W(j);
            currCloneLik2=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),alpha,beta);
            currCloneLik2(aIdx)=bbinopdf_ln(Acounts(aIdx,j),RDmat(aIdx,j),alpha(aIdx),beta(aIdx));
        end
        if i==1 && j==1
            cloneLik1=currCloneLik1;
            cloneLik2=currCloneLik2;
        elseif j==1
            cloneLik1=[cloneLik1 currCloneLik1];
            cloneLik2=[cloneLik2 currCloneLik2];
        else
            cloneLik1(:,i)=cloneLik1(:,i).*currCloneLik1;
            cloneLik2(:,i)=cloneLik2(:,i).*currCloneLik2;
        end
    end
end
%cloneLik1(isnan(expAF))=NaN;
%cloneLik2(isnan(expAF))=NaN;
cloneLik=max(cloneLik1,cloneLik2);
mIdx=cloneLik1>cloneLik2;
pDataSomatic=max(cloneLik,[],2);
cloneId=min(pDataSomatic,0);
alleleId=min(pDataSomatic,0);
for i=1:size(f,2)
    cloneId(pDataSomatic==cloneLik(:,i))=i;
    alleleId(pDataSomatic==cloneLik(:,i) & mIdx(:,i))=1;
    alleleId(pDataSomatic==cloneLik(:,i) & ~mIdx(:,i))=2;
end
%pDataSomatic(Bcounts<1)=min(pDataSomatic(Bcounts<1),pDataHom(Bcounts<1));

for j=1:inputParam.sampleCount
    if inputParam.NormalSample>0
        shiftAF=(cnaF(:,inputParam.NormalSample).*N./M+2*(1-cnaF(:,inputParam.NormalSample))).*cnaF(:,j).*(M./N).*germAF+(1-cnaF(:,j)).*germAF;
        shiftAF(M(:,inputParam.NormalSample)==0,j,1)=0;
        shiftAF=[shiftAF (cnaF(:,inputParam.NormalSample).*N./(N-M)+2*(1-cnaF(:,inputParam.NormalSample))).*cnaF(:,j).*((N-M)./N).*germAF+(1-cnaF(:,j)).*germAF];
        shiftAF(~isfinite(shiftAF(:,1)),1)=germAF(~isfinite(shiftAF(:,1)));
        shiftAF(~isfinite(shiftAF(:,2)),2)=germAF(~isfinite(shiftAF(:,2)));
        currpDataNonDip1=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),min(shiftAF(:,1),1-shiftAF(:,1)).*W(j),max(shiftAF(:,1),1-shiftAF(:,1)).*W(j));
        currpDataNonDip2=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),min(shiftAF(:,2),1-shiftAF(:,2)).*W(j),max(shiftAF(:,2),1-shiftAF(:,2)).*W(j));
    else
        currpDataNonDip1=bbinopdf_ln(Bcounts(:,j),RDmat(:,j),germAF.*W(j),(1-germAF).*W(j));
        currpDataNonDip2=currpDataNonDip1;
    end
    if j==1
        pDataNonDip1=currpDataNonDip1;
        pDataNonDip2=currpDataNonDip2;
    else
        pDataNonDip1=pDataNonDip1.*currpDataNonDip1;
        pDataNonDip2=pDataNonDip2.*currpDataNonDip2;
    end
end

%%%find most likely clone and allele

pDataNonDipMax=max(pDataNonDip1,pDataNonDip2);
pDataNonDipMax(isnan(pDataNonDipMax))=0;

% %%%%find maximum somatic likelihood
% pDataSomatic=NaN(size(Acounts));
% expAFsom=NaN(size(Acounts));
% for i=1:size(f,2)
%     for k=1:2
%         pDataSomatic(cloneId==i & alleleId==k,:)=squeeze(cloneLik(cloneId==i & alleleId==k,i,:,k));
%         expAFsom(cloneId==i & alleleId==k,:)=squeeze(expAF(cloneId==i & alleleId==k,i,:,k));
%     end
% end

%%%% find posterior probabilites
pData=prod(pDataHom,2).*priorHom+prod(pDataHet,2).*priorHet+pDataSomatic.*priorSomatic+prod(pDataOther,2).*priorOther+prod(pDataNonDipMax,2).*priorNonDip;
pSomatic=prod(pDataSomatic,2).*priorSomatic./pData;
pGermline=prod(pDataHet,2).*priorHet./pData;
pHom=prod(pDataHom,2).*priorHom./pData;
pOther=prod(pDataOther,2).*priorOther./pData;
pNonDip=pDataNonDipMax.*priorNonDip./pData;
postComb=array2table([posTable.Chr_1 posTable.Pos_1],'VariableNames',{'Chr', 'Pos'});
postComb.Somatic=pSomatic;
postComb.Het=pGermline;
postComb.Hom=pHom;
postComb.Other=pOther;
postComb.NonDip=pNonDip;
%postComb=array2table([posTable.Chr_1 posTable.Pos_1 pSomatic pGermline pHom pOther pNonDip],'VariableNames',{'Chr', 'Pos', 'Somatic','Het','Hom','Other','NonDip'});
% pDataComb=cell(size(Tcell));
% for i=1:size(Acounts,2)
%     pDataComb{i}=array2table([Tcell{1}.Chr Tcell{1}.Pos pDataSomatic(:,i) pDataHet(:,i) pDataHom(:,i) pDataOther(:,i) pDataNonDipMax(:,i)],'VariableNames',{'Chr', 'Pos', 'Somatic','Het','Hom','Other','NonDip'});
% end
pDataSum=sum(pData);
            
