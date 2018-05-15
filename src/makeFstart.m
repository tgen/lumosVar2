function fStart=makeFstart(fOld,tIdx,inputParam)
%makeFstart - get random starting values for sample fractions
% Syntax:  fStart=makeFstart(fOld,tIdx,inputParam)
%
% Inputs:
%    fOld - sample fractions of current clones
%    tIdx - index of tumor samples
%    inputParam - structure of parameters
%
% Outputs:
%   fStart: matrix of starting values for sample fraction where each column
%       corresponds to a potential clonal variant group and each row corresponds to a sample 
%
% Other m-files required: none
% Other requirements: none
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain, fitCNAmulti

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 7-May-2018

%------------- BEGIN CODE --------------


m=1;
distMat=ones(size(fOld,1),1000);
fStart=[];
    
%%% get matrix of random values of f from a beta distribution with a
%%% mode of priorF for each sample, then find the column of fRand that is
%%% the least similar to the current sample fractions and add it to fStart,
%%% continue until maximum difference is less than minStartDiff
while m/max(distMat(:))>inputParam.minStartDiff
    fRand=100*betarnd(ones(length(tIdx),1000)*inputParam.alphaF,((inputParam.alphaF-1)./inputParam.priorF(tIdx)-inputParam.alphaF+2)*ones(1,1000),size(fOld,1),1000);
    distMat=squareform(pdist([fOld fStart fRand]','cityblock'));
    minList=min(distMat(:,1:size(fOld,2)+size(fStart,2)),[],2);
    [m,idx]=max(minList(size(fOld,2)+size(fStart,2)+1:end));
    fStart=[fStart fRand(:,idx)];
end