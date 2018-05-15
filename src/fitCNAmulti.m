function [segsTable, W, f, CNAscale, nll, t]=fitCNAmulti_v3(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,fInit,cInit,wInit,dbCounts)
%fitCNAmulti - uses EM to fit copy number parameters and estimate copy number
%
% Syntax: [segsTable, W, f, CNAscale, nll, t]=fitCNAmulti(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,fInit,cInit,wInit,dbCounts)
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
%   fInit: initial value for sample fraction matrix
%   wInit: initial value for W parameter
%   dbCounts: counts of common variant positions per segment
%   
% Outputs:
%   segsTable: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio',
%       5-'N',6-'M',7-'F',8-'W',9-'log2FC'
%   W: vector of length of number of samples, controls width of allele
%       frequency distritibutions
%   f: matrix of sample fractions with one row per tumor sample, and one
%       column per clonal variant group
%   CNAscale: centering constant
%   nll: negative log likelihood
%   t: runtimes of optimizations
%
% Other m-files required: nllCNAaddClone, nllCNAmulti, callCNAmulti,
%   joinSNV
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

%%%start parallel pool
pool=gcp('nocreate');
if isempty(pool)
    delete(gcp('nocreate'));
    pc = parcluster('local');
    pc.NumWorkers = inputParam.numCPU;
    parpool(pc, pc.NumWorkers);
end

%%%setup
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
opts2=optimoptions('fmincon','Display','iter','UseParallel',true,'TolX',1e-1,'TolFun',0.1*length(tIdx)./inputParam.addCloneWeight);
j=inputParam.numClones;

%%%%initial optimization
[param{j}, nll(j)]=fmincon(@(param)nllCNAmulti(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,param,dbCounts),[100.*cInit(:); wInit(:); 100*fInit(:)],[],[],[],[],[50*cInit(:); inputParam.minW*ones(size(wInit(:))); zeros(size(fInit(:)));],[200*cInit(:); inputParam.maxW*ones(size(wInit(:))); 100*ones(size(fInit(:)));],[],opts2);
bestMin=nll(j)+j*length(tIdx)./inputParam.addCloneWeight;
prevMin=bestMin;
CNAscale=param{j}(1:length(Tcell))./100
Wcurr=param{j}(length(Tcell)+1:2*length(Tcell))
fOld=reshape(param{j}(2*length(Tcell)+1:end),[],inputParam.numClones)
cBest=CNAscale;
wBest=Wcurr;
fBest=fOld;
iterStuck=0;
while iterStuck<inputParam.iterNoImp
    
    %%%% test recentering
    foundMin=0;
    removeClone=1;
    cTest=CNAscale*ones(1,21);
    cTest(tIdx,:)=(2.^(max(fOld./100,[],2)*[-1:0.1:1])).*(CNAscale(tIdx)*ones(1,21));
    parfor k=1:size(cTest,2)
        nllRecenter(k)=nllCNAmulti(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,[100.*cTest(:,k); Wcurr(:); fOld(:)],dbCounts);
    end
    [~,cIdx]=min(nllRecenter);
    if(any(cTest(:,cIdx)~=CNAscale))
        [paramRecenter, nllRecenter(cIdx)]=fmincon(@(param)nllCNAmulti(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,param,dbCounts),[100.*cTest(:,cIdx); Wcurr(:); fOld(:)],[],[],[],[],[50*cTest(:,cIdx); inputParam.minW*ones(size(wInit(:))); zeros(size(fOld(:)));],[200*cTest(:,cIdx); inputParam.maxW*ones(size(wInit(:))); 100*ones(size(fOld(:)));],[],opts2);
        currMin=nllRecenter(cIdx)+j*length(tIdx)./inputParam.addCloneWeight;
        message=['Recentering']
        CNAscale=paramRecenter(1:length(Tcell))./100
        Wcurr=paramRecenter(length(Tcell)+1:2*length(Tcell))
        fOld=reshape(paramRecenter(2*length(Tcell)+1:end),[],inputParam.numClones)
        if currMin<bestMin
            foundMin=1;
            bestMin=currMin
            cBest=CNAscale;
            wBest=Wcurr;
            fBest=fOld;
        else
            prevMin=currMin;
        end
    end
    
    %%%% test removing clones
    while removeClone && j>=2
        j=size(fOld,2);
        inputParam.numClones=j-1;
        nllRemove=nan(j,1);
        parfor i=1:j
            idx=setdiff(1:j,i);
            fRemove{i}=fOld(:,idx);
            nllRemove(i)=nllCNAmulti(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,[100.*CNAscale(:); Wcurr(:); fRemove{i}(:)],dbCounts);
        end
        [~,removeIdx]=min(nllRemove)
        j=j-1;
        [paramRemove, nll(j)]=fmincon(@(param)nllCNAmulti(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,param,dbCounts),[100.*CNAscale(:); Wcurr(:); fRemove{removeIdx}(:)],[],[],[],[],[50*cInit(:); inputParam.minW*ones(size(wInit(:))); zeros(size(fRemove{removeIdx}(:)));],[200*cInit(:); inputParam.maxW*ones(size(wInit(:))); 100*ones(size(fRemove{removeIdx}(:)));],[],opts2);
        currMin=nll(j)+j*length(tIdx)./inputParam.addCloneWeight;
        if currMin<bestMin
            foundMin=1;
            removeClone=1;
            bestMin=currMin
            CNAscale=paramRemove(1:length(Tcell))./100
            Wcurr=paramRemove(length(Tcell)+1:2*length(Tcell))
            fOld=reshape(paramRemove(2*length(Tcell)+1:end),[],inputParam.numClones)
            cBest=CNAscale;
            wBest=Wcurr;
            fBest=fOld;
        elseif currMin<prevMin
            removeClone=1;
            prevMin=currMin
            CNAscale=paramRemove(1:length(Tcell))./100
            Wcurr=paramRemove(length(Tcell)+1:2*length(Tcell))
            fOld=reshape(paramRemove(2*length(Tcell)+1:end),[],inputParam.numClones)
        else
            j
            message=['not removing any at ' num2str(j)] 
            removeClone=0;
        end
    end  
    
    %%%test adding clones
    addClone=1;
    while addClone
        j=size(fOld,2)+1;
        inputParam.numClones=j;
        pts=makeFstart(fOld,tIdx,inputParam);
        tic
        nllPTS=NaN(size(pts,2),1);
        parfor i=1:size(pts,2)
            nllPTS(i)=nllCNAaddClone(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,CNAscale,Wcurr,fOld,pts(:,i),dbCounts);
        end
        [~,idx]=min(nllPTS)
        toc
        t(j,1)=toc
        fOld
	pts(:,idx)
	fAdd=[fOld pts(:,idx)]
        tic
        [param{j}, nll(j)]=fmincon(@(param)nllCNAmulti(hetPos,somPos,dbPos,Tcell,exonRD,segsMerged,inputParam,param,dbCounts),[100.*CNAscale(:); Wcurr(:); fAdd(:)],[],[],[],[],[50*cInit(:); inputParam.minW*ones(size(wInit(:))); zeros(size(fAdd(:)));],[200*cInit(:); inputParam.maxW*ones(size(wInit(:))); 100*ones(size(fAdd(:)));],[],opts2);
        currMin=nll(j)+j*length(tIdx)./inputParam.addCloneWeight;
        t(j,2)=toc;
        if currMin<bestMin
            foundMin=1;
            addClone=1;
            bestMin=currMin
            CNAscale=param{j}(1:length(Tcell))./100
            Wcurr=param{j}(length(Tcell)+1:2*length(Tcell))
            fOld=reshape(param{j}(2*length(Tcell)+1:end),[],inputParam.numClones)
            cBest=CNAscale;
            wBest=Wcurr;
            fBest=fOld;
        else
            message=['stop adding at ' num2str(j)]
            addClone=0;
            CNAscale=param{j}(1:length(Tcell))./100
            Wcurr=param{j}(length(Tcell)+1:2*length(Tcell))
            fOld=reshape(param{j}(2*length(Tcell)+1:end),[],inputParam.numClones)
            prevMin=currMin
        end
    end 
    if foundMin==0
        iterStuck=iterStuck+1
    end
end

%%%get copy number with optimized parameters
W=wBest
f=fBest./100
CNAscale=cBest
paramOpt=[CNAscale(:); W(:); f(:)];
inputParam.numClones=size(f,2);
[N, M, Ftable, log2FC, cnaIdx]=callCNAmulti(hetPos,Tcell,exonRD,segsMerged,inputParam,paramOpt,dbPos,dbCounts);
segsTable=array2table(segsMerged(:,1:3),'VariableNames',{'Chr','StartPos','EndPos'});
segsTable.N=N;
segsTable.M=M;
segsTable.F=Ftable;
segsTable.cnaIdx=cnaIdx;
for i=1:length(Tcell)
    segsTable.F(:,i)=Ftable(:,i);
    segsTable.log2FC(:,i)=log2FC(:,i);
end

return;
