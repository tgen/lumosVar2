function [segsTable, W, f, CNAscale, nll, t]=fitCNAmulti(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam)
%fitCNA - uses EM to fit copy number parameters and estimate copy number
%
% Syntax: [segsTable, W, f, c, nll, pCNA]=fitCNA(dataHet,dataSom,exonRD,segs,inputParam)
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
%   inputParam: structure with fields: minHetAF, numClones
%   
% Outputs:
%   segsTable: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio',
%       5-'N',6-'M',7-'F',8-'W',9-'log2FC'
%   W: vector of lenght inputParm.numClones, controls width of allele
%       frequency distributions
%   f: vector of sample fraction of each clone
%   c: centering constant
%   nll: negative log likelihood
%   pCNA: vector of probability of copy number state by exon
%
% Other m-files required: nllCNA.m, callCNA.m
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

%%%optimize parameters
%maxClones=inputParam.numClones;
%pool=gcp('nocreate');
%if isempty(pool)
    %delete(gcp('nocreate'));
 %   distcomp.feature( 'LocalUseMpiexec', true);
  %  parpool(inputParam.numCPU);
%end

parfor i=1:length(Tcell)
    diploidPos=(Tcell{i}.BCountF+Tcell{i}.BCountR)./Tcell{i}.ReadDepthPass>inputParam.minHetAF;
    cInit(i,:)=median(2*Tcell{i}.ControlRD(diploidPos & hetPos)./Tcell{i}.ReadDepthPass(diploidPos & hetPos));
    wInit(i,:)=nanmedian(exonRD{i}(:,4));
    dataHet=[Tcell{i}.Chr(hetPos) Tcell{i}.Pos(hetPos) Tcell{i}.ControlRD(hetPos) Tcell{i}.ReadDepthPass(hetPos) Tcell{i}.BCountF(hetPos)+Tcell{i}.BCountR(hetPos)];
    if i==inputParam.NormalSample
        dataSom=[];
    else
        dataSom=[Tcell{i}.Chr(somPos) Tcell{i}.Pos(somPos) Tcell{i}.ControlRD(somPos) Tcell{i}.ReadDepthPass(somPos) Tcell{i}.BCountF(somPos)+Tcell{i}.BCountR(somPos)];
    end
    paramInit=fminsearchbnd(@(param)nllCNA(dataHet,dataSom,exonRD{i},segsMerged,inputParam,param),[cInit(i,:) wInit(i,:) 0.5],[cInit(i,:)*0.5 inputParam.minW 0],[cInit(i,:)*2 inputParam.maxW 1]);
    CNAscale(i,:)=paramInit(1)
    Wcurr(i,:)=paramInit(2)
    fInit(i,:)=100.*paramInit(3)
end

tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
CNAscale
Wcurr
fInit=fInit(tIdx,:)
fOld=[];
%chi2p(1)=0;
opts=optimoptions('fmincon','TolX',1e-1,'TolFun',1e-1,'Display','none');
opts2=optimoptions('fmincon','Display','iter','UseParallel',true,'TolX',1e-2,'TolFun',1e-2);
j=1;
%while max(chi2p)<0.05
while 1
    inputParam.numClones=j;
  %ms=MultiStart('Display','iter','UseParallel',true);
    if j==1
        pts=[fInit makeFstart(fInit,inputParam)];
    else
        pts=makeFstart(fOld,inputParam);
    end
    %problem=createOptimProblem('fmincon','objective',@(fNew)nllCNAaddClone(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam,cInit,Wcurr,fOld,fNew),'x0',pts(:,1),'lb',zeros(size(fInit)),'ub',100*ones(size(fInit)),'options',opts);
    %startPoints=CustomStartPointSet(pts');
    %paramMS=run(ms,problem,startPoints);
    tic
    parfor i=1:size(pts,2)
        if isfinite(nllCNAaddClone(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam,cInit,wInit,fOld,pts(:,i)))
            [paramPTS{i}, nllPTS(i)]=fmincon(@(fNew)nllCNAaddClone(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam,cInit,wInit,fOld,fNew),pts(:,i),[],[],[],[],zeros(size(fInit)),100*ones(size(fInit)),[],opts);
        else
            message=['not defined at ' num2str(pts(:,i))]
            nllPTS(i)=INF;
        end
    end
    [minNLL,idx]=min(nllPTS)
    t(j,1)=toc
    fOld=[fOld paramPTS{idx}]
    tic
    [param{j}, nll(j)]=fmincon(@(param)nllCNAmulti_v2(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam,param),[100.*cInit(:); wInit(:); fOld(:)],[],[],[],[],[50*cInit(:); inputParam.minW*ones(size(wInit(:))); zeros(size(fOld(:)));],[200*cInit(:); inputParam.maxW*ones(size(wInit(:))); 100*ones(size(fOld(:)));],[],opts2);
    t(j,2)=toc;
    CNAscale=param{j}(1:length(Tcell))./100
    Wcurr=param{j}(length(Tcell)+1:2*length(Tcell))
    fOld=reshape(param{j}(2*length(Tcell)+1:end),[],inputParam.numClones)
    if j>1
        %    chi2p(j)=1-chi2cdf(2*(nll(j-1)-nll(j)),length(Tcell)-1)
        if inputParam.addCloneWeight*(nll(j-1)-nll(j))<length(tIdx)
            inputParam.numClones=j-1;
            for i=1:size(fOld,2)
                idx=setdiff(1:j,i);
                fRemove=fOld(:,idx);
                [paramRemove{i}, nllRemove(i)]=fmincon(@(param)nllCNAmulti_v2(hetPos,somPos,Tcell,exonRD,segsMerged,inputParam,param),[100.*cInit(:); wInit(:); fRemove(:)],[],[],[],[],[50*cInit(:); inputParam.minW*ones(size(wInit(:))); zeros(size(fRemove(:)));],[200*cInit(:); inputParam.maxW*ones(size(wInit(:))); 100*ones(size(fRemove(:)));],[],opts2);
            end
            [removeNLL,removeIdx]=min(nllRemove)
            CNAscale=paramRemove{removeIdx}(1:length(Tcell))./100
            Wcurr=paramRemove{removeIdx}(length(Tcell)+1:2*length(Tcell))
            fOld=reshape(paramRemove{removeIdx}(2*length(Tcell)+1:end),[],inputParam.numClones)
            if nll(j-1)<=removeNLL
                break;
            else
                nll(j-1)=removeNLL;
                continue;
            end
        end
    end
    j=j+1;
end

%%%use optimized parameters to call copy number
%cloneIdx=find(chi2p<0.05,1,'last');
%cloneIdx=j-1;

CNAscale=paramRemove{removeIdx}(1:length(Tcell))./100
W=Wcurr
f=fOld./100
paramOpt=[CNAscale(:); W(:); f(:)];
%inputParam.numClones=cloneIdx;
[N, M, Ftable, log2FC, cnaIdx]=callCNAmulti(hetPos,Tcell,exonRD,segsMerged,inputParam,paramOpt);
segsTable=array2table(segsMerged(:,1:3),'VariableNames',{'Chr','StartPos','EndPos'});
segsTable.N=N;
segsTable.M=M;
segsTable.F=Ftable;
segsTable.cnaIdx=cnaIdx;

for i=1:length(Tcell)
    segsTable.F(:,i)=Ftable(:,i);
    %segsTable.W(:,i)=Wtable(:,i);
    segsTable.log2FC(:,i)=log2FC(:,i);
end

return;
