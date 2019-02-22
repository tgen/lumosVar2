function writeCloneSummary(segsTable,exonRD,Tcell,fIn,cloneId,inputParam,Filter,somaticDetected,sampleFrac)
%writeCloneSummary - writes summary tables and plots
%
% Syntax:  writeCloneSummary(segsTable,exonRD,Tcell,fIn,cloneId,inputParam,Filter,somaticDetected,sampleFrac)
%
% Inputs:
%   segsTable: table of segment data with variables: {'Chr','StartPos','EndPos',
%       'N','M','F','cnaIdx'}
%   exonRD: cell array of matrices of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   Tcell: cell array of tables with columns: {'Chr','Pos','ReadDepthPass'}
%   fIn: sample fraction matrix
%   cloneId: index of clonal variant group
%   inputParam: structure of parameters
%   Filter: cell array of variant calls
%   somaticDetected: logical matrix indicating which samples somatic
%       variants were detected in
%   sampleFrac: estimated fraction of cells containing each variant
%
% Outputs:
%    output file names/paths specified by inputParam.outName
%    *.cloneSummary.tsv - table clonal variant group data
%    *.cloneSummary.pdf - plots of clonal variant groups
%    *.groupLinePlots.pdf - plots of sampleFractions of somatic variants in
%    each group
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
% Last revision: 9-May-2018
%------------- BEGIN CODE --------------

%%%fill in sample fractions for normal sample
if inputParam.NormalSample>0
    f=[zeros(inputParam.sampleCount,inputParam.numClones) ones(inputParam.sampleCount,1)];
    tIdx=setdiff(1:inputParam.sampleCount,inputParam.NormalSample);
    f(tIdx,1:end-1)=[fIn];
else
    f=[fIn ones(inputParam.sampleCount,1)];
end

%%%create table of variant counts by clonal group
sampleNames=char(regexp(inputParam.sampleNames,',','split')');
sampleNamesShort=cellstr(sampleNames(:,1:min(namelengthmax-3,size(sampleNames,2))));
cloneTable=array2table([f'; NaN(1,size(f,1))],'VariableNames',strcat('SF_',regexprep(cellstr(sampleNamesShort),'-|\.','_')));
cloneTable{size(f,2),1:size(f,1)}=1-max(f(:,1:end-1),[],2)';
cloneTable.Row=[cellstr(strcat('ClonalGroup_',num2str([1:size(f,2)-1]'))); 'Germline'; 'Pair'];
T=Tcell{1};
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

%%% find exon cloneId by exon
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

%%% write summary table
writetable([cloneTable array2table([CNcount; NaN(1,size(CNcount,2))],'VariableNames',CNname)],[inputParam.outName '.cloneSummary.tsv'],'Delimiter','\t','WriteRowNames',1,'FileType','text');

%%% summary plot of sample fractions (y) by sample (x) where line thickness
%%% is proportional to the number of copy number alterations in clonal
%%% variant groups and size of circle is proportional to the number of
%%% small mutations in the clonal variant group
colors=linspecer(size(f,2)-1);
subplot(3,2,1);
hold on;
for i=1:size(colors,1)
    plot(f(:,i),'-','color',colors(i,:),'LineWidth',10*(sum(CNcount(i,:))+1)./(sum(CNcount(:))+1));
    scatter([1:size(f,1)],f(:,i),400*(cloneTable.somaticDetected(i,:)+1)./sum(cloneTable.somaticPass+1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',colors(i,:));
end
ylim([0 1]);
set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
ylabel('Clonal Sample Fraction');

%%% if there is more than one sample, lineplot of sample fraction vs
%%% sample, where each line is a somatic variant and color indicates clonal
%%% variant group
if (size(sampleFrac,2)>1)
    subplot(3,2,3);
    hold on;
    for i=1:size(colors,1)
        plot(sampleFrac(strcmp(Filter,'SomaticPASS') & cloneId(:,1)==i,:)','Color',colors(i,:));
    end
    ylim([0 1]);
    set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
    ylabel('Variant Sample Fraction');
end

%%%beeswarm plot of sample fractions by sample
if(size(sampleFrac,2)>1)
    subplot(3,2,5);
else
    subplot(3,2,3);
end
cloneIdplot=cloneId(:,1);
cloneIdplot(~strcmp(Filter,'SomaticPASS'))=NaN;
cloneIdplot(strcmp(Filter,'SomaticPairPASS'))=0;
[~,gName]=grp2idx(cloneIdplot);
cIdx=str2double(gName);
cIdx=cIdx(cIdx>0);
somIdx=~isnan(cloneIdplot);
if sum(strcmp(Filter,'SomaticPairPASS'))<1
    plotSpread(min(sampleFrac(somIdx,:),1),'CategoryIdx',repmat(cloneIdplot(somIdx),1,size(sampleFrac,2)),'CategoryColors',colors(cIdx,:));
else
    plotSpread(min(sampleFrac(somIdx,:),1),'CategoryIdx',repmat(cloneIdplot(somIdx),1,size(sampleFrac,2)),'CategoryColors',[0 0 0; colors(cIdx,:)])
end
set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
ylabel('Variant Sample Fraction');

%%%barplot of number of copy number altered exons by copy number state
subplot(3,2,2)
b=bar(CNcount','stacked');
for i=1:size(colors,1)
    b(i).FaceColor=colors(i,:);
end
b(end).FaceColor=[0 0 0];
set(gca,'XTick',1:length(CNname),'XTickLabel',CNname,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',90);
axis tight;
legend([cellstr(num2str([1:size(colors,1)]')); {'Germ'}]);
xlabel('Copy Number State');
ylabel('Number of Exons');

%%%barplot of number of somatic variants detected by sample
subplot(3,2,4)
b=bar([cloneTable.somaticDetected'; zeros(1,height(cloneTable))],'stacked');
for i=1:size(colors,1)
    b(i).FaceColor=colors(i,:);
end
if inputParam.NormalSample>0
    b(end).FaceColor=[0 0 0];
end
set(gca,'XTick',1:size(sampleNames,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
axis tight;
xl=get(gca,'xlim');
xlim([xl(1) xl(2)-1]);
ylabel('Number of Somatic Variants Detected');

%%%if more than one sample, barplot of the counts of common and unique
%%%variants detected in combinations of samples
if(size(sampleFrac,2)>1 & sum(strcmp(Filter,'SomaticPASS')|strcmp(Filter,'SomaticPairPASS'))>0)
    subplot(3,2,6);
    detectCat=cellstr(num2str(somaticDetected));
    [gIdx,gName]=grp2idx(detectCat);
    cIdx=strcmp(Filter,'SomaticPASS');
    h2=histcounts2(gIdx(cIdx),cloneId(cIdx,1),1:length(gName)+1,1:max(cloneId(:,1))+1);
    h=histcounts(gIdx(strcmp(Filter,'SomaticPairPASS')),1:length(gName)+1);
    hcomb=[h2 h'];
    [~,ord]=sort(sum(hcomb,2),'descend');
    b=bar(hcomb(ord,:),'stacked');
    for i=1:size(h2,2)
        b(i).FaceColor=colors(i,:);
    end
    b(end).FaceColor=[0 0 0];
    gNum=cell2mat(cellfun(@(x) str2double(strsplit(x)),gName,'UniformOutput',0));
    axis tight;
    ylm=ylim();
    hold on;
    for i=1:size(gNum,2)
        x=find(gNum(ord,i));
        scatter(x,-0.1*ylm(2)*i*ones(size(x)),'.k');
    end
    xlm=[0.5 sum(sum(hcomb(ord,:),2)>0)+0.5];
    xlim(xlm);
    text(ones(size(sampleNamesShort'))*xlm(2),linspace(-0.1*ylm(2),-0.1*ylm(2)*length(sampleNamesShort),length(sampleNamesShort)),sampleNamesShort,'HorizontalAlignment','left','Interpreter','none','FontSize',8);
    yt=get(gca,'YTick');
    set(gca,'YTick',yt(yt>0),'FontSize',8,'XTick',[]);
    legend([cellstr(num2str([1:size(colors,1)]')); {'NA'}]);
end
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 7.5 10]);
print(gcf,'-dpdf',[inputParam.outName '.cloneSummary.pdf'],'-r300');
close(gcf);

%%% if more than one sample, seperate line plots of each clonal variant
%%% group where lines are colored by mean read depth
if(size(sampleFrac,2)>1)
    RD=zeros(height(Tcell{1}),1);
    for i=1:inputParam.sampleCount
        RD=RD+Tcell{i}.ReadDepthPass;
    end
    RD=RD./inputParam.sampleCount;
    RDbin=ceil(log2(RD+1));
    cmap=jet(max(RDbin));
    
    for i=1:size(colors,1)
        subplot(ceil((size(colors,1)+2)/2),2,i);
        currIdx=strcmp(Filter,'SomaticPASS') & cloneId(:,1)==i;
        hold on;
         for j=1:max(RDbin)
            plot(min(sampleFrac(currIdx & RDbin==j,:)',1),'-','Color',cmap(j,:));
        end
        ylim([0 1]);
        set(gca,'XTick',1:size(f,1),'XTickLabel',sampleNames,'FontSize',8,'TickLabelInterpreter','none','XTickLabelRotation',10);
        ylabel('Variant Sample Fraction');
        title(['Variant Group ' num2str(i)]);
    end
    subplot(ceil((size(colors,1)+2)/2),2,i+1);
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
