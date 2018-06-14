function nll = nllCNAmulti(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,param,dbCounts)
%nllCNAmulti - finds the negative log likelihood of the data given the
%               parameters
%
% Syntax: nll = nllCNAmulti(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,param,dbCounts)
%
% Inputs:
%   hetPos: logical vector indicating heterozygous positions
%   somPos: logical vector indicating somatic positions
%   dbPos: logical vector indicating common germline variant positions
%   Tcell - cell array of tables with length equal to the number of bams,
%       each table must have the following columns:  {'Chr','Pos', 'ReadDepthPass',
%       'RefComb','AComb','AcountsComb','AmeanBQ','BComb','BCountsComb','BmeanBQ',
%       'ApopAFcomb','BpopAFcomb','CosmicCount','PosMapQC'};      
%   exonRD: cell array oi matrices of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segsMerged: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%   inputParam: structure with fields: minHetAF, numClones
%   param: 1:numSamples - CNAscale
%          numSamples+1:2*numSamples - W
%          2*numSamples+1:end - f
%   dbCounts: counts of common variant positions per segment
%   
% Outputs:
%   nll - negative log likelihood of model given parameters
%
% Other m-files required: callCNAmulti, joinSNV
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 8-May-2018

%------------- BEGIN CODE --------------

%%% read in inputs
E=table();
E.Chr=exonRD{1}(:,1);
E.StartPos=exonRD{1}(:,2);
E.EndPos=exonRD{1}(:,3);
for i=1:length(Tcell)
    E.TumorRD(:,i)=exonRD{i}(:,4);
    E.NormalRD(:,i)=exonRD{i}(:,5);
end
CNAscale=param(1:length(Tcell))./100;
W=param(length(Tcell)+1:2*length(Tcell));
f=reshape(param(2*length(Tcell)+1:end)./100,[],inputParam.numClones);

%%%get copy number calls
[N, M, Ftable, ~, cnaIdx, nllCNA]=callCNAmulti(hetPos,Tcell,exonRD,segsMerged,inputParam,[CNAscale(:); W(:); f(:)],dbPos,dbCounts);
segsTable=array2table(segsMerged(:,1:3),'VariableNames',{'Chr','StartPos','EndPos'});
segsTable.N=N;
segsTable.M=M;
segsTable.F=Ftable;
segsTable.cnaIdx=cnaIdx;

%%%get likelihoods of somatic variants
Tsom=cell(size(Tcell));
for j=1:length(Tcell)
    T=Tcell{j}(somPos,:);
    idx=getPosInRegionSplit([T.Chr T.Pos],segsTable{:,1:3},inputParam.blockSize);
    T.NumCopies=segsTable.N(idx);
    T.MinAlCopies=segsTable.M(idx);
    T.cnaF=segsTable.F(idx,j);
    Tsom{j}=T;
end
somLik=NaN(height(Tsom{1}),length(Tsom));
somIdx=NaN(height(Tsom{1}),length(Tsom));
[~, ~,pDataComb,clones,~]=jointSNV(Tsom, f, W, inputParam);
for j=1:length(Tsom)
    somLik(:,j)=pDataComb{j}.Somatic+inputParam.minLik;
    somIdx(:,j)=clones;
end

%%%find likelihood of association between copy number state and clones
cnState=strcat(strsplit(sprintf('%d\n',Tsom{1}.NumCopies)),'_',strsplit(sprintf('%d\n',Tsom{1}.MinAlCopies)));
cnState=cnState(1:end-1);
if (height(Tsom{1})>0)
   [~,~,chiP,~]=crosstab(somIdx(:,1),cnState);
   chiP=min(chiP,1);
else
   chiP=1;
end

totalPosCount=sum(exonRD{1}(:,3)-exonRD{1}(:,2));

%%%find clonal fraction difference priors
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
priorF=ones(sum(somPos),1);
if(sum(somPos)>0)
    fDiff=nan(size(f,2),1);
    for i=1:size(f,2)
        if size(f,1)>1
            fDiff(i)=geomean(pdist(f(:,i),'cityblock')+0.05);
        else
            fDiff(i)=geomean(pdist([0; f(:,i); 1],'cityblock')+0.05);
        end
    end
    if length(inputParam.priorF)>1
        priorFDiff=geomean(pdist(inputParam.priorF,'cityblock')+0.05);
    else
        priorFDiff=geomean(pdist([0; inputParam.priorF; 1],'cityblock')+0.05);
    end
    for j=1:length(fDiff)
        priorF(somIdx(:,1)==j)=betapdf(fDiff(j),inputParam.alphaF,(inputParam.alphaF-1)./priorFDiff-inputParam.alphaF+2)+inputParam.minLik;
    end
else
    somLik=1;
    priorF=1;
end

%%%find negative log likliehood
nll=-(sum((log(priorF)+sum(log(somLik),2)))./size(somLik,1)+nllCNA+log(chiP+realmin)./(inputParam.priorSomaticSNV*totalPosCount));

return;
