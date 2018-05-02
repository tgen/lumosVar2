function [Filter,passPos,somaticDetected,trustScore,artifactScore]=callVariants(Tcell,P,inputParam)

if inputParam.NormalSample>0
    bIdx=Tcell{inputParam.NormalSample}.AcountsComb>=Tcell{inputParam.NormalSample}.BcountsComb;
    aIdx=Tcell{inputParam.NormalSample}.AcountsComb<Tcell{inputParam.NormalSample}.BcountsComb;
else
    bIdx=Tcell{1}.ApopAFcomb>=Tcell{1}.BpopAFcomb;
    aIdx=Tcell{1}.ApopAFcomb<Tcell{1}.BpopAFcomb;
end

somaticDetected=zeros(size(Tcell{1},1),length(Tcell));
altCount=zeros(size(Tcell{1},1),length(Tcell));
for i=1:length(Tcell)
    T=Tcell{i};
   if inputParam.NormalSample<1
        somaticDetected(P.Somatic(:,i)>0.5 & P.DataSomatic(:,i)>P.DataHom(:,i) & (T.BcountsComb>0 | T.A~=T.RefComb),i)=1;
   else
        somaticDetected(P.SomaticPair(:,i)>0.5 | (P.Somatic(:,i)>0.5 & P.DataSomatic(:,i)>P.DataHom(:,i) & (T.BcountsComb>0 | T.A~=T.RefComb)),i)=1;
   end
   altCount(somaticDetected(:,i) & bIdx,i)=T.BcountsComb(somaticDetected(:,i) & bIdx);
   altCount(somaticDetected(:,i) & aIdx,i)=T.AcountsComb(somaticDetected(:,i) & aIdx);
   altCount(max([P.Somatic P.SomaticPair],[],2)<0.5 & P.Het(:,i)>P.Hom(:,i),i)=T.BcountsComb(max([P.Somatic P.SomaticPair],[],2)<0.5 & P.Het(:,i)>P.Hom(:,i));
   altCount(max([P.Somatic P.SomaticPair],[],2)<0.5 & P.Hom(:,i)>P.Het(:,i),i)=T.AcountsComb(max([P.Somatic P.SomaticPair],[],2)<0.5 & P.Hom(:,i)>P.Het(:,i));
end
%[mCount,mIdx]=max(altCount,[],2);
trustScore=1-prod((1-P.trust+realmin).^altCount,2).^(1./sum(altCount,2));
%trustScore=sum((altCount+1).*P.trust,2)./sum(altCount+1,2);
%artifactScore=sum((altCount+1).*P.artifact,2)./sum(altCount+1,2);
artifactScore=prod((P.artifact+realmin).^altCount,2).^(1./sum(altCount,2));
passPos=trustScore>inputParam.pGoodThresh & artifactScore<inputParam.pGoodThresh;

%%% assign filters
% passPos=max(P.trust,[],2)>inputParam.pGoodThresh & max(P.artifact,[],2)<inputParam.pGoodThresh;
% for i=1:length(Tcell)
%     pos=mCount>0 & mIdx==i;
%     passPos(pos)=P.trust(pos,i)>inputParam.pGoodThresh & P.artifact(pos,i)<inputParam.pGoodThresh;
% end
%trustSom=P.trust;
%trustSom(~somaticDetected)=0;
Filter=cell(size(passPos));
Filter(P.Somatic(:,1)>0.5,:)={'SomaticLowQC'};
Filter(P.Somatic(:,1)>inputParam.pSomaticThresh & passPos,:)={'SomaticPASS'};
%Filter(P.Somatic(:,1)>0.5 & min([Tcell{1}.ApopAFcomb Tcell{1}.BpopAFcomb],[],2)>inputParam.maxSomPopFreq,:)={'SomaticDBsnp'};
Filter(P.Het(:,1)>0.5,:)={'GermlineHetLowQC'};
Filter(P.Het(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHetPASS'};
Filter(P.Hom(:,1)>0.5,:)={'GermlineHomLowQC'};
Filter(P.Hom(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHomPASS'};
Filter(P.NonDip(:,1)>0.5,:)={'GermlineShiftLowQC'};
Filter(P.NonDip(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineShiftPASS'};
Filter(P.Somatic(:,1)<0.5 & P.Het(:,1)<0.5 & P.Hom(:,1)<0.5 & P.NonDip(:,1)<0.5,:)={'NoCall'};
%Filter(min(P.artifact,[],2)>inputParam.pGoodThresh,:)={'REJECT'};
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
if inputParam.NormalSample>0
    idxSom=strncmp(Filter,'Somatic',7);
    %idxSomPass=strncmp(Filter,'SomaticPASS',11);
    for i=1:length(tIdx)
        %Filter(P.SomaticPair(:,tIdx(i))>0.5 & ~idxSom & min(P.artifact(:,[tIdx(i) inputParam.NormalSample]),[],2)<inputParam.pGoodThresh,:)={'SomaticPairLowQC'};
        Filter(P.SomaticPair(:,tIdx(i))>0.5 & ~idxSom,:)={'SomaticPairLowQC'};
        Filter(P.SomaticPair(:,tIdx(i))>inputParam.pSomaticThresh & ~idxSom & passPos)={'SomaticPairPASS'};
    end
end
idxSom=strncmp(Filter,'Somatic',7);
Filter(idxSom & passPos & min([Tcell{1}.ApopAFcomb Tcell{1}.BpopAFcomb],[],2)>inputParam.maxSomPopFreq,:)={'SomaticDBsnp'};
Filter(artifactScore>inputParam.pGoodThresh,:)={'REJECT'};
