function nll = nllCNAmulti_v4(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,param,dbCounts)
%nllCNA - computes negative loglikliehood of copy number parameters
%
% Syntax: nll = nllCNA(dataHet,dataSom,exonRD,segs,inputParam,param)
%
% Inputs:
%   dataHet: data for germline heterozygous positions with columns:
%       1-'Chr',2-'Pos',3-'ControlRD',4-'TumorRD',5-'Bcount'
%   dataSom: data for somatic positions with columns:
%       1-'Chr',2-'Pos',3-'ControlRD',4-'TumorRD',5-'Bcount'
%   exonRD: matrix of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segs: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%   inputParam: structure with fields: cnaPrior, minAllelePrior, minLik,
%       alphaF, priorF
%   param: vector of paramaters of length 2*numClones+1
%       param(1) - copy number scaling constant
%       param(2:numClones+1)=W (controls width of allele frequency dist)
%       param(numClones+1:2*numClones+1)=f (sample fractions)
%   
% Outputs:
%   nll: negative log likelihood
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, fitCNA, callCNA

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

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
    
[N, M, Ftable, ~, cnaIdx, nllCNA]=callCNAmulti_v3(hetPos,Tcell,exonRD,segsMerged,inputParam,[CNAscale(:); W(:); f(:)],dbPos,dbCounts);
segsTable=array2table(segsMerged(:,1:3),'VariableNames',{'Chr','StartPos','EndPos'});
segsTable.N=N;
segsTable.M=M;
segsTable.F=Ftable;
segsTable.cnaIdx=cnaIdx;

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
[~, ~,pDataComb,clones,~]=jointSNV_v2(Tsom, f, W, inputParam);
for j=1:length(Tsom)
    %[~,locb]=ismember([Tsom{j}.Chr Tsom{j}.Pos],[postComb.Chr postComb.Pos],'rows');
    somLik(:,j)=pDataComb{j}.Somatic+inputParam.minLik;
    somIdx(:,j)=clones;
end

cnState=strcat(strsplit(sprintf('%d\n',Tsom{1}.NumCopies)),'_',strsplit(sprintf('%d\n',Tsom{1}.MinAlCopies)));
cnState=cnState(1:end-1);

if (height(Tsom{1})>0)
   [~,~,chiP,~]=crosstab(somIdx(:,1),cnState);
   chiP=min(chiP,1);
else
   chiP=1;
end

% 
% for j=1:length(Tcell)
%     T=Tcell{j}(dbPos,:);
%     idx=getPosInRegionSplit([T.Chr T.Pos],segsTable{:,1:3},inputParam.blockSize);
%     T.NumCopies=segsTable.N(idx);
%     T.MinAlCopies=segsTable.M(idx);
%     T.cnaF=segsTable.F(idx,j);
%     Tdb{j}=T;
% end
% postComb=jointSNV_v2(Tdb, f, W, inputParam);
% somDBpos=postComb.Somatic>inputParam.pSomaticThresh;
totalPosCount=sum(exonRD{1}(:,3)-exonRD{1}(:,2));
% [~,p]=fishertest([sum(somDBpos) inputParam.dbSNPposCount-sum(somDBpos); sum(somPos) totalPosCount-sum(somPos)],'tail','right');



%%% find likelihood of somatic variant
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
priorF=ones(sum(somPos),1);
%%% find likelihood of somatic variant
if(sum(somPos)>0)
    [fMax,fIdx]=max(f,[],1);
    %for j=1:length(tIdx)
        %priorF(:,tIdx(j))=betapdf(f(j,somIdx(:,tIdx(j))),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF(tIdx(j))-inputParam.alphaF+2)+inputParam.minLik;
            %end
    for j=1:length(fMax)
        priorF(somIdx(:,1)==j)=betapdf(fMax(j),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF(tIdx(fIdx(j)))-inputParam.alphaF+2)+inputParam.minLik;
    end
else
    somLik=1;
    priorF=1;
end
    
%end
%priorF=1;

%%% sum negative log likliehoods
%nll=sum((-sum(log(somLik))-sum(log(hetlikMax))-sum(log(depthlikMax))-sum(log(priorCNAMax))-sum(log(priorMinAlleleMax))-sum(log(priorF))-nansum(log(priorCNAf)))./(length(somLik)+length(hetlikMax)+length(depthlikMax)+length(priorCNAMax)+length(priorMinAlleleMax)+length(priorF)+sum(~isnan(priorCNAf))));
%nll=sum((-sum(log(somLik))-sum(log(hetlikMax))-sum(log(depthlikMax))))./(length(somLik)+length(hetlikMax)+length(depthlikMax));

%nll=-sum((log(priorF)+sum(log(somLik),2))./(inputParam.priorSomaticSNV*sum(E.EndPos-E.StartPos)))+nllCNA-log(p+realmin)./((sum(somPos)./totalPosCount)*inputParam.dbSNPposCount)+log(chiP+realmin)./(inputParam.priorSomaticSNV*totalPosCount);
%sum((log(priorF)+sum(log(somLik),2)))./(inputParam.priorSomaticSNV*totalPosCount)
%sum((log(priorF)+sum(log(somLik),2)))./size(somLik,2);
%nllCNA;
%log(p+realmin)./(inputParam.priorSomaticSNV*inputParam.dbSNPposCount)
%log(chiP+realmin)./(inputParam.priorSomaticSNV*totalPosCount);
%nll=-(sum((log(priorF)+sum(log(somLik),2)))./(inputParam.priorSomaticSNV*totalPosCount)+nllCNA+log(p+realmin)./(inputParam.priorSomaticSNV*inputParam.dbSNPposCount)+log(chiP+realmin)./(inputParam.priorSomaticSNV*totalPosCount));
%nll=-(sum((log(priorF)+sum(log(somLik),2)))./(inputParam.priorSomaticSNV*totalPosCount)+nllCNA+log(chiP+realmin)./(inputParam.priorSomaticSNV*totalPosCount));
nll=-(sum((log(priorF)+sum(log(somLik),2)))./size(somLik,1)+nllCNA+log(chiP+realmin)./(inputParam.priorSomaticSNV*totalPosCount));

return;
