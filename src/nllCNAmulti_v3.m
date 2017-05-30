function nll = nllCNAmulti_v3(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam,param)
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
    
[N, M, Ftable, ~, cnaIdx, nllCNA]=callCNAmulti(hetPos,Tcell,exonRD,segsMerged,inputParam,[CNAscale(:); W(:); f(:)]);
segsTable=array2table(segsMerged(:,1:3),'VariableNames',{'Chr','StartPos','EndPos'});
segsTable.N=N;
segsTable.M=M;
segsTable.F=Ftable;
segsTable.cnaIdx=cnaIdx;

for j=1:length(Tcell)
    T=Tcell{j}(somPos,:);
    idx=getPosInRegionSplit([T.Chr T.Pos],segsTable{:,1:3},inputParam.blockSize);
    T.NumCopies=segsTable.N(idx);
    T.MinAlCopies=segsTable.M(idx);
    T.cnaF=segsTable.F(idx,j);
    Tsom{j}=T;
end

[postComb, ~,pDataComb,clones,~]=jointSNV_v2(Tsom, f, W, inputParam);
for j=1:length(Tsom)
    [lia,locb]=ismember([Tsom{j}.Chr Tsom{j}.Pos],[postComb.Chr postComb.Pos],'rows');
    somLik(:,j)=pDataComb{j}.Somatic(locb);
    somIdx(:,j)=clones(locb);
end


%%% find likelihood of somatic variant
if(sum(somPos)>0)
    fMax=max(f,[],1);
    priorF=betapdf(fMax(somIdx),inputParam.alphaF,(inputParam.alphaF-1)./inputParam.priorF-inputParam.alphaF+2)+inputParam.minLik;
else
    somLik=1;
    priorF=1;
end
    
%end
%priorF=1;

%%% sum negative log likliehoods
%nll=sum((-sum(log(somLik))-sum(log(hetlikMax))-sum(log(depthlikMax))-sum(log(priorCNAMax))-sum(log(priorMinAlleleMax))-sum(log(priorF))-nansum(log(priorCNAf)))./(length(somLik)+length(hetlikMax)+length(depthlikMax)+length(priorCNAMax)+length(priorMinAlleleMax)+length(priorF)+sum(~isnan(priorCNAf))));
%nll=sum((-sum(log(somLik))-sum(log(hetlikMax))-sum(log(depthlikMax))))./(length(somLik)+length(hetlikMax)+length(depthlikMax));

nll=sum(-sum(log(somLik)./(inputParam.priorSomaticSNV*sum(E.EndPos-E.StartPos)))-mean(log(priorF)))+nllCNA;

return;