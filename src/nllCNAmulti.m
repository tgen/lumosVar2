function nll = nllCNAmulti(hetData,somData,dbData,exonRD,segsMerged,inputParam,param,dbCounts)
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
for i=1:inputParam.sampleCount
    E.TumorRD(:,i)=exonRD{i}(:,4);
    E.NormalRD(:,i)=exonRD{i}(:,5);
end
CNAscale=param(1:inputParam.sampleCount)./100;
W=param(inputParam.sampleCount+1:2*inputParam.sampleCount);
f=reshape(param(2*inputParam.sampleCount+1:end)./100,[],inputParam.numClones);

%%%get copy number calls
[N, M, Ftable, ~, cnaIdx, nllCNA]=callCNAmulti(hetData,exonRD,segsMerged,inputParam,[CNAscale(:); W(:); f(:)],dbData,dbCounts);
segsTable=array2table(segsMerged(:,1:3),'VariableNames',{'Chr','StartPos','EndPos'});
segsTable.N=N;
segsTable.M=M;
segsTable.F=Ftable;
segsTable.cnaIdx=cnaIdx;

%%%get likelihoods of somatic variants
somData.N=segsTable.N(somData.idxSeg);
somData.M=segsTable.M(somData.idxSeg);
for j=1:inputParam.sampleCount
    somData.(strcat('cnaF_',num2str(i-1)))=segsTable.F(somData.idxSeg,j);
end
somLik=NaN(height(somData),inputParam.sampleCount);
somIdx=NaN(height(somData),inputParam.sampleCount);
[~, ~, somIdx,~,~,somLik]=jointSNV(somData, f, W, inputParam);

%%%find likelihood of association between copy number state and clones
cnState=strcat(strsplit(sprintf('%d\n',somData.N)),'_',strsplit(sprintf('%d\n',somData.M)));
cnState=cnState(1:end-1);
if (height(somData)>0)
   [~,~,chiP,~]=crosstab(somIdx(:,1),cnState);
   chiP=min(chiP,1);
else
   chiP=1;
end

totalPosCount=sum(exonRD{1}(:,3)-exonRD{1}(:,2));

%%%find clonal fraction difference priors
tIdx=setdiff(1:inputParam.sampleCount,inputParam.NormalSample);
priorF=ones(height(somData),1);
if(height(somData)>0)
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
