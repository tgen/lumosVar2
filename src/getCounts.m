function [countsAll]=getCounts(Tcell, inputParam)
%getCounts - finds two alleles with highest counts across all samples 
% Syntax:  countsAll=getCounts(paramFile)
%
% Inputs:
%    Tcell - cell array of tables with length equal to the number of bams,
%       each table must have the following columns:  {'Chr','Pos','ReadDepth',
%       'ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ','AmeanMQ',
%       'AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ','BmeanMQ',
%       'BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount','ControlRD',
%        'PosMapQC','perReadPass','abFrac'};
%    inputParam - data structure with fields from paramFile   
%
% Outputs:
%   countsAll - table with the following variables:
%     'Chr' - chromosome
%     'Pos' - positition on chr
%     'Ref' - reference allele
%     'A' - allele the most counts across all samples
%     'B' - allele the second most counts across all samples
%     'ApopAF' - population allele frequency of A allele 
%     'BpopAF' - population allele frequency of B allele 
%     'cosmicCount' - maximum cosmic count observed at position
%     'Acounts' - matrix with one column per sample with read counts for A
%     'Bcounts' - matrix with one column per sample with read counts for A
%
%
% Other m-files required: none
% Other requirements: none
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 4-May-2018

%------------- BEGIN CODE --------------

%%%find unique positions
posList=[];
for i=1:size(Tcell,2)
    T=Tcell{i};
    posList=[posList; T.Chr T.Pos];
end
posList=gather(unique(posList,'rows'));

%%%read data from Tcell into matrices
numVar=size(posList,1);
%Amat=NaN(numVar,size(Tcell,2));
%Bmat=NaN(numVar,size(Tcell,2));
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
    BmeanBQ(locb,i)=T.BmeanBQ;
    AmeanBQ(locb,i)=T.AmeanBQ;
    mapQC(locb,i)=T.PosMapQC;
end

%%%%set values when A and B are the same across samples
Ref=zeros(size(Amat));
Ref(min(RefOrig,[],2)==max(RefOrig,[],2),:)=RefOrig(min(RefOrig,[],2)==max(RefOrig,[],2),:);
for i=1:size(BmeanBQ,2)
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

%%% find values for A and B when they are different across samples
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

countsAll=array2table([posList Ref A B ApopAF BpopAF cosmic],'VariableNames',{'Chr','Pos','Ref','A','B','ApopAF','BpopAF','cosmicCount'});
for i=1:size(Acounts,2)
    countsAll.Acounts(:,i)=Acounts(:,i);
    countsAll.Bcounts(:,i)=Bcounts(:,i);
end
    

            
