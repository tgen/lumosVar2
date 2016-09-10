function writeCloneSummary(segsTable,Ecell,Tcell,P,fIn,cloneId,inputParam)
%writeCloneSummary - writes summary of variant counts by clone
%
% Syntax: writeCloneSummary(segsTable,E,T,pSomatic,posterior,f,W,cloneId,inputParam)
%
% Inputs:
%   segsTable: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio',
%       5-'N',6-'M',7-'F',8-'W',9-'log2FC'
%    E - table of exon data with the following columns: 'Chr','StartPos','EndPos',
%       'TumorRD','NormalRD', 'MapQC', 'perReadPass', 'abFrac'
%    T - table of position data with the following columns: 'Chr','Pos',
%       'ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ',
%       'AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ',
%       'BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount',
%       'ControlRD','PosMapQC','perReadPass','abFrac'
%   pSomatic: posterior probability somatic variant
%   pTrust: posterior probability call should be trusted
%   pArtifact: posterior probability artificat
%   f: sample fraction of each clone
%   W: w parameter for each clone
%   cloneId: most likley clone assuming somatic
%   inputParam: structure with all parameters   
%
% Outputs:
%    writes a csv file
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, callCNA, callSNV

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016
%------------- BEGIN CODE --------------

f=zeros(length(Tcell),inputParam.numClones);
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
f(tIdx,1:end)=reshape(fIn,[],inputParam.numClones);

cloneTable=array2table(f','VariableNames',regexp(inputParam.sampleNames,',','split'));

T=Tcell{1};
rejectPos=max(P.artifact,[],2)>inputParam.pGoodThresh;
passPos=max(P.trust,[],2)>inputParam.pGoodThresh & T.RefComb>0 & T.Acomb>0 & T.Bcomb>0 & ~rejectPos;
%%% count variants by clone
for i=1:size(cloneTable,1)
    pass(i,:)=sum(cloneId(:,1)==i & passPos & P.Somatic(:,1)>inputParam.pSomaticThresh & min([T.ApopAF T.BpopAF],[],2)<inputParam.maxSomPopFreq);
    lowqc(i,:)=sum(cloneId(:,1)==i & ~rejectPos & P.Somatic(:,1)>0.5 & min([T.ApopAF T.BpopAF],[],2)<inputParam.maxSomPopFreq)-pass(i,:);
    db(i,:)=sum(cloneId(:,1)==i & ~rejectPos & P.Somatic(:,1)>0.5 & min([T.ApopAF T.BpopAF],[],2)>=inputParam.maxSomPopFreq);
end
cloneTable.somaticPass=pass;
cloneTable.somaticLowQC=lowqc;
cloneTable.somaticDB=db;

message=['made clone table']

%%% find exon cloneId
E=Ecell{1};
idx=getPosInRegions([E.Chr mean([E.StartPos E.EndPos],2)],segsTable{:,1:3});
exonCN(~isnan(idx),:)=[segsTable.N(idx(~isnan(idx))) segsTable.M(idx(~isnan(idx))) segsTable.F(idx(~isnan(idx)),1)];
exonCloneId=zeros(length(exonCN),1);
for i=1:size(f,2)
    exonCloneId(exonCN(:,3)==f(1,i) & (exonCN(:,1)~=2 | exonCN(:,2)~=1))=i;
end

message=['found exon clone id']

%%% count CNV by clone
n=1;
for j=0:4
    for k=0:floor(j/2)
        CNname(n)={['N' num2str(j) '_M' num2str(k)]};
        for i=1:size(f,2)
            CNcount(i,n)=sum(exonCloneId==i & exonCN(:,1)==j & exonCN(:,2)==k);
        end
        n=n+1;
    end
end
CNname(n)={'N5toInf'};
for i=1:size(f,2)
    CNcount(i,n)=sum(exonCloneId==i & exonCN(:,1)>=5);
end
message=['about summary table']

writetable([cloneTable array2table(CNcount,'VariableNames',CNname)],[inputParam.outName '.cloneSummary.csv']);

cmap=colormap('hsv');
colors=cmap(1:64./size(f,2):64,:);
hold on;
for i=1:size(f,2)
    plot(f(:,i),'-o','MarkerSize',50*(cloneTable.somaticPass(i)+1)./sum(cloneTable.somaticPass),'Color',colors(i,:),'MarkerFaceColor',colors(i,:),'LineWidth',10*(sum(CNcount(i,:))+1)./sum(CNcount(:)));
end
set(gca,'XTick',1:size(f,2),'XTickLabel',cloneTable.Properties.VariableNames(1:size(f,2)),'FontSize',14);
ylabel('Variant Sample Fraction');

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 7 4]);
print(gcf,'-dpng',[inputParam.outName '.cloneSummary.png'],'-r300');
close(gcf);
