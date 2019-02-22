function [Filter,passPos,somaticDetected,trustScore]=callVariants(posTable,P,inputParam)
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
aIdx=posTable.aIdx;
bIdx=~posTable.aIdx;


%%%%count alternate alleles observed for each variant in each sample
somaticDetected=P.Somatic>0.5 & ~P.homPos(:,1) | P.SomaticPair(:,1)>0.5;
altCount=posTable.BcountsComb_0;
altCount(aIdx & P.Hom<0.5 | P.Hom>=0.5 & ~aIdx)=posTable.AcountsComb_0(aIdx & P.Hom<0.5 | P.Hom>=0.5 & ~aIdx);
for i=2:inputParam.sampleCount
   somaticDetected=[somaticDetected (P.Somatic>0.5 & ~P.homPos(:,i) | P.SomaticPair(:,i)>0.5)];
   altCount=[altCount posTable.(strcat('BcountsComb_',num2str(i-1)))];
   altCount(aIdx & P.Hom<0.5 | P.Hom>=0.5 & ~aIdx,i)=posTable.(strcat('AcountsComb_',num2str(i-1)))(aIdx & P.Hom<0.5 | P.Hom>=0.5 & ~aIdx);
end
%%%calculate trust and artifact scores weighting by alt counts
trustScore=1-prod((1-min(P.trust,1)+realmin).^altCount,2).^(1./sum(altCount,2));
%artifactScore=prod((max(P.artifact,0)+realmin).^altCount,2).^(1./sum(altCount,2));
passPos=trustScore>inputParam.pGoodThresh;

%%%determine variant calls
filterNames={'SomaticLowQC','SomaticPASS','GermlineHetLowQC','GermlineHetPASS','GermlineHomLowQC','GermlineHomPASS','GermlineShiftLowQC','GermlineShiftPASS','NoCall','SomaticPairLowQC','SomaticPairPASS','SomaticDBsnp','REJECT'};
Filter=min(posTable.Chr_0,0);
Filter(P.Somatic(:,1)>0.5,:)=1;
Filter(P.Somatic(:,1)>inputParam.pSomaticThresh & passPos,:)=2;
Filter(P.Het(:,1)>0.5,:)=3;
Filter(P.Het(:,1)>inputParam.pGermlineThresh & passPos,:)=4;
Filter(P.Hom(:,1)>0.5,:)=5;
Filter(P.Hom(:,1)>inputParam.pGermlineThresh & passPos,:)=6;
Filter(P.NonDip(:,1)>0.5,:)=7;
Filter(P.NonDip(:,1)>inputParam.pGermlineThresh & passPos,:)=8;
Filter(P.Somatic(:,1)<0.5 & P.Het(:,1)<0.5 & P.Hom(:,1)<0.5 & P.NonDip(:,1)<0.5,:)=9;
idxSom=(Filter==1 | Filter==2); 
if inputParam.NormalSample>0
    for i=1:inputParam.sampleCount
        Filter(P.SomaticPair(:,i)>0.5 & ~idxSom,:)=10;
        Filter(P.SomaticPair(:,i)>inputParam.pSomaticThresh & ~idxSom & passPos)=11;
    end
end
Filter(idxSom & passPos & min([posTable.ApopAF_0 posTable.BpopAF_0],[],2)>inputParam.maxSomPopFreq,:)=12;
Filter(1-trustScore>inputParam.pGoodThresh,:)=13;
%Filter(artifactScore>inputParam.pGoodThresh,:)=13;
Filter=categorical(Filter,[1:13],filterNames);
