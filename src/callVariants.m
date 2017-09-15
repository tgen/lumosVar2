function [Filter,passPos,somaticDetected]=callVariants(Tcell,P,inputParam)

somaticDetected=zeros(size(Tcell{1},1),length(Tcell));
for i=1:length(Tcell)
    T=Tcell{i};
   if inputParam.NormalSample<1
        somaticDetected(P.Somatic(:,i)>0.5 & P.DataSomatic(:,i)>P.DataHom(:,i) & (T.BcountsComb>0 | T.A~=T.RefComb),i)=1;
   else
        somaticDetected(P.SomaticPair(:,i)>0.5 | (P.Somatic(:,i)>0.5 & P.DataSomatic(:,i)>P.DataHom(:,i) & (T.BcountsComb>0 | T.A~=T.RefComb)),i)=1;
   end
end

%%% assign filters
passPos=max(P.trust,[],2)>inputParam.pGoodThresh & max(P.artifact,[],2)<inputParam.pGoodThresh & Tcell{1}.RefComb>0 & Tcell{1}.Acomb>0 & Tcell{1}.Bcomb>0;
trustSom=P.trust;
trustSom(~somaticDetected)=0;
Filter=cell(size(passPos));
Filter(P.Somatic(:,1)>0.5,:)={'SomaticLowQC'};
Filter(P.Somatic(:,1)>inputParam.pSomaticThresh & max(trustSom,[],2)>inputParam.pGoodThresh & passPos,:)={'SomaticPASS'};
Filter(P.Somatic(:,1)>0.5 & min([Tcell{1}.ApopAF Tcell{1}.BpopAF],[],2)>inputParam.maxSomPopFreq,:)={'SomaticDBsnp'};
Filter(P.Het(:,1)>0.5,:)={'GermlineHetLowQC'};
Filter(P.Het(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHetPASS'};
Filter(P.Hom(:,1)>0.5,:)={'GermlineHomLowQC'};
Filter(P.Hom(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHomPASS'};
Filter(P.NonDip(:,1)>0.5,:)={'GermlineShiftLowQC'};
Filter(P.NonDip(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineShiftPASS'};
Filter(P.Somatic(:,1)<0.5 & P.Het(:,1)<0.5 & P.Hom(:,1)<0.5 & P.NonDip(:,1)<0.5,:)={'NoCall'};
Filter(min(P.artifact,[],2)>inputParam.pGoodThresh,:)={'REJECT'};
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample)
if inputParam.NormalSample>0
    idxSom=strncmp(Filter,'Somatic',7);
    idxSomPass=strncmp(Filter,'SomaticPASS',11);
    for i=1:length(tIdx)
        Filter(P.SomaticPair(:,tIdx(i))>0.5 & ~idxSom & min(P.artifact(:,[tIdx(i) inputParam.NormalSample]),[],2)<inputParam.pGoodThresh,:)={'SomaticPairLowQC'};
        Filter(max(P.SomaticPair(:,tIdx(i)),[],2)>inputParam.pSomaticThresh & ~idxSomPass & P.trust(:,tIdx(i))>inputParam.pGoodThresh &  max(P.artifact(:,[tIdx(i) inputParam.NormalSample]),[],2)<inputParam.pGoodThresh,:)={'SomaticPairPASS'};
    end
end
