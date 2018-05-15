function [Filter,passPos,somaticDetected,trustScore,artifactScore]=callVariants(Tcell,P,inputParam)
% takes posterior probabilities of variant classes and variant qualities
% and classifies variants 
%
% Syntax: [Filter,passPos,somaticDetected,trustScore,artifactScore]=callVariants(Tcell,P,inputParam)
% 
% Inputs:
%   Tcell: cell arry of tables each with the following columns: 'Chr','Pos',
%       'AcountsComb', 'BcountsComb', 'ApopAFcomb', 'BpopAFcomb'
%   P: table of probabilites with the following variables: 'Somatic',
%   'SomaticPair', 'Het', 'Hom', 'trust', 'artifact', 'DataSomatic',
%   'DataHom', 'NonDip'
%   
% Outputs:
%    Filter - variant calls
%    passPos - logical vector of pass positions
%    somaticDetected - matrix of positions where somatic variants are
%    detected
%    trustScore - score indicating if position is trusted
%    artifactScore -score indicating if position is artifact
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain, jointSNV, qualDiscrimCalls

%%%determine whether A or B allele is the candidate variant
if inputParam.NormalSample>0
    bIdx=Tcell{inputParam.NormalSample}.AcountsComb>=Tcell{inputParam.NormalSample}.BcountsComb;
    aIdx=Tcell{inputParam.NormalSample}.AcountsComb<Tcell{inputParam.NormalSample}.BcountsComb;
else
    bIdx=Tcell{1}.ApopAFcomb>=Tcell{1}.BpopAFcomb;
    aIdx=Tcell{1}.ApopAFcomb<Tcell{1}.BpopAFcomb;
end

%%%%count alternate alleles observed for each variant in each sample
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
%%%calculate trust and artifact scores weighting by alt counts
trustScore=1-prod((1-P.trust+realmin).^altCount,2).^(1./sum(altCount,2));
artifactScore=prod((P.artifact+realmin).^altCount,2).^(1./sum(altCount,2));
passPos=trustScore>inputParam.pGoodThresh & artifactScore<inputParam.pGoodThresh;

%%%determine variant calls
Filter=cell(size(passPos));
Filter(P.Somatic(:,1)>0.5,:)={'SomaticLowQC'};
Filter(P.Somatic(:,1)>inputParam.pSomaticThresh & passPos,:)={'SomaticPASS'};
Filter(P.Het(:,1)>0.5,:)={'GermlineHetLowQC'};
Filter(P.Het(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHetPASS'};
Filter(P.Hom(:,1)>0.5,:)={'GermlineHomLowQC'};
Filter(P.Hom(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHomPASS'};
Filter(P.NonDip(:,1)>0.5,:)={'GermlineShiftLowQC'};
Filter(P.NonDip(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineShiftPASS'};
Filter(P.Somatic(:,1)<0.5 & P.Het(:,1)<0.5 & P.Hom(:,1)<0.5 & P.NonDip(:,1)<0.5,:)={'NoCall'};
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
if inputParam.NormalSample>0
    idxSom=strncmp(Filter,'Somatic',7);
    for i=1:length(tIdx)
        Filter(P.SomaticPair(:,tIdx(i))>0.5 & ~idxSom,:)={'SomaticPairLowQC'};
        Filter(P.SomaticPair(:,tIdx(i))>inputParam.pSomaticThresh & ~idxSom & passPos)={'SomaticPairPASS'};
    end
end
idxSom=strncmp(Filter,'Somatic',7);
Filter(idxSom & passPos & min([Tcell{1}.ApopAFcomb Tcell{1}.BpopAFcomb],[],2)>inputParam.maxSomPopFreq,:)={'SomaticDBsnp'};
Filter(artifactScore>inputParam.pGoodThresh,:)={'REJECT'};
