function writeCloneSummary(segsTable,exonRD,Tcell,fIn,cloneId,inputParam,Filter,somaticDetected)
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
sampleNamesShort=cellstr(sampleNames(:,1:min(namelengthmax-2,size(sampleNames,2))));
cloneTable=array2table([f'; NaN(1,size(f,1))],'VariableNames',strcat('S_',regexprep(cellstr(sampleNamesShort),'-','_')));

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
%E=Ecell{1};
idx=getPosInRegions([exonRD{1}(:,1) mean(exonRD{1}(:,2:3),2)],segsTable{:,1:3});
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
%if inputParam.NormalSample>0
colors=linspecer(size(f,2)-1);
%else
%    colors=linspecer(size(f,2));
%end
subplot(3,2,1);
hold on;
for i=1:size(colors,1)
    plot(f(:,i),'-','color',colors(i,:),'LineWidth',10*(sum(CNcount(i,:))+1)./sum(CNcount(:)));
    %if rem(i,2)==0
    %plot(f(:,i),'-o','MarkerSize',50*(cloneTable.somaticDetected(i,:)+1)./sum(cloneTable.somaticPass+1),'color',colors(i,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',colors(i,:),'LineWidth',10*(sum(CNcount(i,:))+1)./sum(CNcount(:)));
    scatter([1:size(f,1)],f(:,i),400*(cloneTable.somaticDetected(i,:)+1)./sum(cloneTable.somaticPass+1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',colors(i,:));
    %else
    %plot(f(:,i),'-o','MarkerSize',50*(cloneTable.somaticDetected(i,:)+1)./sum(cloneTable.somaticPass+1),'color',colors(i,:),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',colors(i,:),'LineWidth',10*(sum(CNcount(i,:))+1)./sum(CNcount(:)));
    %   scatter([1:size(f,1)],f(:,i),400*(cloneTable.somaticDetected(i,:)+1)./sum(cloneTable.somaticPass+1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',colors(i,:));
    %end
end
ylim([0 1]);
set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
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

if (size(sampleFrac,2)>1)
    subplot(3,2,3);
    hold on;
    for i=1:size(colors,1)
        plot(sampleFrac(strcmp(Filter,'SomaticPASS') & cloneId(:,1)==i,:)','Color',colors(i,:));
        %plot(sampleFrac(strcmp(Filter,'SomaticPairPASS') & cloneId(:,1)==i,:)',':','Color',colors(i,:));
    end
    %plot(sampleFrac(strcmp(Filter,'SomaticPairPASS'),:)','Color','k');
    ylim([0 1]);
    set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
    ylabel('Variant Sample Fraction');
end

if(size(sampleFrac,2)>1)
    subplot(3,2,5);
else
    subplot(3,2,3);
end
% plot(sampleFrac(strcmp(Filter,'SomaticPairPASS'),:)');
% ylim([0 1]);
% set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
% ylabel('Variant Sample Fraction');
cloneIdplot=cloneId(:,1);
cloneIdplot(~strcmp(Filter,'SomaticPASS'))=NaN;
cloneIdplot(strcmp(Filter,'SomaticPairPASS'))=0;
[~,gName]=grp2idx(cloneIdplot);
%gName
%colors
cIdx=str2double(gName);
cIdx=cIdx(cIdx>0);
somIdx=~isnan(cloneIdplot);
if inputParam.NormalSample<1
    plotSpread(min(sampleFrac(somIdx,:),1),'CategoryIdx',repmat(cloneIdplot(somIdx),1,size(sampleFrac,2)),'CategoryColors',colors);
else
    plotSpread(min(sampleFrac(somIdx,:),1),'CategoryIdx',repmat(cloneIdplot(somIdx),1,size(sampleFrac,2)),'CategoryColors',[0 0 0; colors(cIdx,:)])
end
set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
ylabel('Variant Sample Fraction');

subplot(3,2,2)
b=bar(CNcount','stacked');
for i=1:size(colors,1)
    b(i).FaceColor=colors(i,:);
end
%if inputParam.NormalSample>0
b(end).FaceColor=[0 0 0];
%end

set(gca,'XTick',1:length(CNname),'XTickLabel',CNname,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',90);
axis tight;
legend([cellstr(num2str([1:size(colors,1)]')); {'Germ'}]);
xlabel('Copy Number State');
ylabel('Number of Exons');

subplot(3,2,4)
b=bar([cloneTable.somaticDetected'; zeros(1,height(cloneTable))],'stacked');
%b=bar(cloneTable.somaticDetected','stacked');
for i=1:size(colors,1);
    b(i).FaceColor=colors(i,:);
end
if inputParam.NormalSample>0
    b(end).FaceColor=[0 0 0];
end
set(gca,'XTick',1:size(sampleNames,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
%xlim([0.5 length(sampleNames)+0.5]);
axis tight;
xl=get(gca,'xlim');
xlim([xl(1) xl(2)-1]);
%legend(num2str([1:size(CNcount,1)]'));
ylabel('Number of Somatic Variants Detected');

if(size(sampleFrac,2)>1)
    subplot(3,2,6);
    detectCat=cellstr(num2str(somaticDetected));
    [gIdx,gName]=grp2idx(detectCat);
    cIdx=strcmp(Filter,'SomaticPASS');
    %[tbl,~,~,labels]=crosstab(gIdx(cIdx),cloneId(cIdx,1));
    h2=histcounts2(gIdx(cIdx),cloneId(cIdx,1),1:length(gName)+1,1:max(cloneId(:,1))+1);
    h=histcounts(gIdx(strcmp(Filter,'SomaticPairPASS')),1:length(gName)+1);
    hcomb=[h2 h'];
    gPos=find(sum(hcomb')>0);
    [~,ord]=sort(sum(hcomb(gPos,:)'),'descend');
    b=bar(hcomb(gPos(ord),:),'stacked');
    for i=1:size(h2,2)
        b(i).FaceColor=colors(i,:);
    end
    b(end).FaceColor=[0 0 0];
    gNum=cell2mat(cellfun(@(x) str2double(strsplit(x)),gName,'UniformOutput',0));
    axis tight;
    ylm=ylim();
    hold on;
    for i=1:size(gNum,2)
        x=find(gNum(gPos(ord),i));
        scatter(x,-0.1*ylm(2)*i*ones(size(x)),'.k');
    end
    xlm=xlim();
    text(ones(size(sampleNamesShort'))*xlm(1),linspace(-0.1*ylm(2),-0.1*ylm(2)*length(sampleNamesShort),length(sampleNamesShort)),sampleNamesShort,'HorizontalAlignment','right','Interpreter','none','FontSize',8);
    yt=get(gca,'YTick');
    set(gca,'YTick',yt(yt>0),'FontSize',8,'XTick',[]);
    legend([cellstr(num2str([1:size(colors,1)]')); {'NA'}]);
end

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 7.5 10]);
print(gcf,'-dpdf',[inputParam.outName '.cloneSummary.pdf'],'-r300');
close(gcf);

if(size(sampleFrac,2)>1)
    RD=zeros(height(Tcell{1}),1);
    for i=1:length(Tcell)
        RD=RD+Tcell{i}.ReadDepthPass;
    end
    RD=RD./length(Tcell);
    RDbin=ceil(log2(RD+1));
    cmap=jet(max(RDbin));
    
    for i=1:size(colors,1)
        subplot(ceil((size(colors,1)+2)/2),2,i);
        currIdx=strcmp(Filter,'SomaticPASS') & cloneId(:,1)==i;
        hold on;
        %scatter([1:length(Tcell)]'*ones(1,sum(currIdx)),sampleFrac(currIdx,:)');
        for j=1:max(RDbin)
            plot(min(sampleFrac(currIdx & RDbin==j,:)',1),'-','Color',cmap(j,:));
        end
        ylim([0 1]);
        set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
        ylabel('Variant Sample Fraction');
        title(['Variant Group ' num2str(i)]);
    end
    subplot(ceil((size(colors,1)+2)/2),2,i+1);
    %plot(sampleFrac(strcmp(Filter,'SomaticPairPASS'),:)','Color','k');
    currIdx=strcmp(Filter,'SomaticPairPASS');
    hold on;
    for j=1:max(RDbin)
        plot(min(sampleFrac(currIdx & RDbin==j,:)',1),'-','Color',cmap(j,:));
    end
    ylim([0 1]);
    set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
    ylabel('Variant Sample Fraction');
    title('Ungrouped Variants');
    subplot(ceil((size(colors,1)+2)/2),2,i+2);
    colormap('jet');
    caxis([0 max(RDbin)]);
    colorbar('northoutside','XTick',0:2:max(RDbin),'XTickLabel',2.^[0:2:max(RDbin)]);
    title('Mean Read Depth');
    axis off;
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [1 1 7.5 10]);
    print(gcf,'-dpdf',[inputParam.outName '.groupLinePlots.pdf'],'-r300');
    close(gcf);
end