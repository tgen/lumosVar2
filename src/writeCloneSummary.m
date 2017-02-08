function writeCloneSummary(segsTable,Ecell,Tcell,fIn,cloneId,inputParam,Filter,somaticDetected)
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
if inputParam.NormalSample>0
    f=[zeros(length(Tcell),inputParam.numClones) ones(length(Tcell),1)];
    tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
    f(tIdx,1:end-1)=[fIn];
else
    f=[fIn ones(length(Tcell),1)];
end

sampleNames=char(regexp(inputParam.sampleNames,',','split')');
%sampleNamesShort=cellstr(sampleNames(:,1:25));
cloneTable=array2table([f'; NaN(1,size(f,1))],'VariableNames',cellstr(sampleNames));

T=Tcell{1};
%rejectPos=max(P.artifact,[],2)>inputParam.pGoodThresh;
%passPos=max(P.trust,[],2)>inputParam.pGoodThresh & T.RefComb>0 & T.Acomb>0 & T.Bcomb>0 & ~rejectPos;
%%% count variants by clone
for i=1:size(f,2)
    pass(i,:)=sum(cloneId(:,1)==i & strcmp(Filter,'SomaticPASS'));
    detect(i,:)=sum(somaticDetected(cloneId(:,1)==i & strcmp(Filter,'SomaticPASS'),:),1);
    lowqc(i,:)=sum(cloneId(:,1)==i & strcmp(Filter,'SomaticLowQC'));
    db(i,:)=sum(cloneId(:,1)==i & strcmp(Filter,'SomaticDBsnp'));
end
cloneTable.somaticPass=[pass; sum(strcmp(Filter,'SomaticPairPASS'))];
cloneTable.somaticDetected=[detect; sum(somaticDetected(strcmp(Filter,'SomaticPairPASS'),:),1)];
cloneTable.somaticLowQC=[lowqc; sum(strcmp(Filter,'SomaticPairLowQC'))];
cloneTable.somaticDB=[db; NaN];


message=['made clone table']

%%% find exon cloneId
E=Ecell{1};
idx=getPosInRegions([E.Chr mean([E.StartPos E.EndPos],2)],segsTable{:,1:3});
exonCN(~isnan(idx),:)=[segsTable.N(idx(~isnan(idx))) segsTable.M(idx(~isnan(idx))) segsTable.F(idx(~isnan(idx)),1)];
exonCloneId(~isnan(idx),:)=segsTable.cnaIdx(idx(~isnan(idx)));
exonCloneId(exonCN(:,1)==2 & exonCN(:,2)==1,:)=0;


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
message=['made summary table']

writetable([cloneTable array2table([CNcount; NaN(1,size(CNcount,2))],'VariableNames',CNname)],[inputParam.outName '.cloneSummary.csv']);

colors=linspecer(size(f,2));
subplot(2,1,1);
hold on;
for i=1:size(f,2)
    if rem(i,2)==0
        plot(f(:,i),'-o','MarkerSize',50*(cloneTable.somaticPass(i)+1)./sum(cloneTable.somaticPass+1),'color',colors(i,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',colors(i,:),'LineWidth',10*(sum(CNcount(i,:))+1)./sum(CNcount(:)));
    else
        plot(f(:,i),'-o','MarkerSize',50*(cloneTable.somaticPass(i)+1)./sum(cloneTable.somaticPass+1),'color',colors(i,:),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',colors(i,:),'LineWidth',10*(sum(CNcount(i,:))+1)./sum(CNcount(:)));
    end
end
ylim([0 1]);
set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',5);
ylabel('Clonal Sample Fraction');


if inputParam.NormalSample>0
    bIdx=Tcell{inputParam.NormalSample}.AcountsComb>=Tcell{inputParam.NormalSample}.BcountsComb;
else
    bIdx=T.ApopAFcomb>=T.BpopAFcomb;
end

for j=1:length(Tcell)
    T=Tcell{j};
    AF(bIdx,j)=T.BcountsComb(bIdx)./T.ReadDepthPass(bIdx);
    AF(~bIdx,j)=T.AcountsComb(~bIdx)./T.ReadDepthPass(~bIdx);
    matchIdx=f(j,cloneId(:,j))'==T.cnaF;
    sampleFrac(matchIdx,j)=AF(matchIdx,j).*(T.cnaF(matchIdx).*T.NumCopies(matchIdx)+2.*(1-T.cnaF(matchIdx)))./(T.NumCopies(matchIdx)-T.MinAlCopies(matchIdx));
    sampleFrac(~matchIdx,j)=AF(~matchIdx,j).*(T.cnaF(~matchIdx).*T.NumCopies(~matchIdx)+2.*(1-T.cnaF(~matchIdx)));
end

subplot(2,1,2);
hold on;
for i=1:size(f,2)
    plot(sampleFrac(strcmp(Filter,'SomaticPASS') & cloneId(:,1)==i,:)','Color',colors(i,:));
    plot(sampleFrac(strcmp(Filter,'SomaticPairPASS') & cloneId(:,1)==i,:)',':','Color',colors(i,:));
end
%plot(sampleFrac(strcmp(Filter,'SomaticPairPASS'),:)','Color','k');
ylim([0 1]);
set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',5);
ylabel('Variant Sample Fraction');
% 
% subplot(3,1,3);
% hold on;
% for i=1:size(f,2)
%     plot(AF(strcmp(Filter,'SomaticPASS') & cloneId(:,1)==i,:)','Color',colors(i,:));
%     plot(AF(strcmp(Filter,'SomaticPairPASS') & cloneId(:,1)==i,:)',':','Color',colors(i,:));
% end
% %plot(sampleFrac(strcmp(Filter,'SomaticPairPASS'),:)','Color','k');
% set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',5);
% ylabel('Variant Sample Fraction');
% 

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 7 4]);
print(gcf,'-dpng',[inputParam.outName '.cloneSummary.png'],'-r300');
close(gcf);
