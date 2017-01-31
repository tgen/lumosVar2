function plotTumorOnly(exonRD,segsTable,CNAscale,fIn,Tcell,somPos,hetPos,cloneId,inputParam)
%plotTumorOnly - plots summary figure of tumor only caller results
%
% Syntax: writeCloneSummary(segsTable,E,T,pSomatic,posterior,f,W,cloneId,inputParam)
%
% Inputs:
%   exonRD: matrix of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segsTable: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio',
%       5-'N',6-'M',7-'F',8-'W',9-'log2FC'
%   f: sample fraction of each clone
%   T - table of position data with the following columns: 'Chr','Pos',
%       'ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ',
%       'AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ',
%       'BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount',
%   somPos: logical index of postions in T that are somatic
%   hetPos: logical index of positions in T that are germline heterozygous
%   cloneId: most likley clone assuming somatic
%   inputParam: structure with all parameters   
%
% Outputs:
%    writes a png file
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
    f=[fIn];
end
%%% transform chromosome coord to linear coord
maxChr=max(segsTable.Chr);
T=Tcell{1};
for i=1:maxChr
    chrLen(i)=max(segsTable{segsTable.Chr==i,3})+1E7;
end
chrOffset=[0; cumsum(chrLen(1:21))'];
for i=1:maxChr
    segCoord(segsTable.Chr==i,1)=(segsTable{segsTable.Chr==i,2}+chrOffset(i))/1E6;
    segCoord(segsTable.Chr==i,2)=(segsTable{segsTable.Chr==i,3}+chrOffset(i))/1E6;
    exonCoord(exonRD{1}(:,1)==i,1)=(exonRD{1}(exonRD{1}(:,1)==i,2)+chrOffset(i))/1E6;
    exonCoord(exonRD{1}(:,1)==i,2)=(exonRD{1}(exonRD{1}(:,1)==i,3)+chrOffset(i))/1E6;
    Tcoord(T.Chr==i,1)=(T.Pos(T.Chr==i)+chrOffset(i))/1E6;
end

sexChr=regexp(inputParam.sexChr,',','split');
chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr']

%%% calculate log2FC
for i=1:length(Tcell)
    log2FC(:,i)=log2((CNAscale(i)./2).*exonRD{i}(:,4)./exonRD{i}(:,5));
    Nlog2R(:,i)=log2(segsTable.F(:,i).*segsTable.N/2+(1-segsTable.F(:,i)));
end

%%% plot exon log2FC
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1]);
    for i=1:2:maxChr
        idx=exonRD{j}(:,1)==i;
        scatter(mean(exonCoord(idx,:),2),log2FC(idx,j),1,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
        hold on;
    end
    for i=2:2:maxChr
        idx=exonRD{j}(:,1)==i;
        scatter(mean(exonCoord(idx,:),2),log2FC(idx,j),1,'.','MarkerFaceColor',[0.8235    0.7059    0.5490],'MarkerEdgeColor',[0.8235    0.7059    0.5490])
        hold on;
    end
    ylabel('log2(FoldChange)','FontSize',10);
end

%cmap=colormap('hsv');
%colors=cmap(1:64./size(f,2):64,:);
colors=linspecer(size(f,2));
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,j*5);
    hold on;
    for i=1:size(f,2)
        b=bar(i,f(j,i));
        b.FaceColor=colors(i,:);
        if(rem(i,2)==0)
            b.EdgeColor=[0 0 0];
        else
            b.EdgeColor=[1 1 1];
        end
    end
    ylim([0 1]);
    xlim([0 size(f,2)+1]);
    ylabel('Sample Fraction','FontSize',10);
end

%%% plot segments
samples=regexp(inputParam.sampleNames,',','split');
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1])
    for i=1:size(f,2)
        pos=segsTable.cnaIdx==i & (segsTable.N~=2 | segsTable.M~=1);
        if(rem(i,2)==0)
            plot(segCoord(pos,:)',ones(2,1)*Nlog2R(pos,j)','Color',[0 0 0],'linewidth',4);
            plot(segCoord(pos,:)',ones(2,1)*Nlog2R(pos,j)','Color',colors(i,:),'linewidth',3);
        else
            plot(segCoord(pos,:)',ones(2,1)*Nlog2R(pos,j)','Color',colors(i,:),'linewidth',4); 
        end
        %plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)',colors(i),'linewidth',4);
    end
    pos=(segsTable.N==2 & segsTable.M==1);
    plot(segCoord(pos,:)',ones(2,1)*Nlog2R(pos,j)','k','linewidth',3);
    %plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)','k','linewidth',4);
    %ticks=[0 2.^[1:7]];
    %tickpos=log2(median(segsTable(:,7)).*ticks/2+(1-median(segsTable(:,7))));
    %set(gca,'YTick',tickpos,'YTickLabel',ticks,'tickDir','out');
    set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrList,'FontSize',6);
    axis([0 max(segCoord(:,2)) min(Nlog2R(:))-1 max(Nlog2R(:))+1])
    title(samples{j},'FontSize',10,'Interpreter','none');
end
subplot(length(Tcell)+1,5,[5*(length(Tcell)+1)-4 5*(length(Tcell)+1)-1])
hold on;
for i=1:size(f,2)
    pos=segsTable.cnaIdx(:,1)==i & (segsTable.N~=2 | segsTable.M~=1);
    if(rem(i,2)==0)
        plot(segCoord(pos,:)',log2(ones(2,1)*segsTable.N(pos)'+1),'Color',[0 0 0],'linewidth',3);
        plot(segCoord(pos,:)',log2(ones(2,1)*segsTable.N(pos)'+1),'Color',colors(i,:),'linewidth',2);
    else
        plot(segCoord(pos,:)',log2(ones(2,1)*segsTable.N(pos)'+1),'Color',colors(i,:),'linewidth',3);
    end
    %plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)',colors(i),'linewidth',4);
end
pos=(segsTable.N==2 & segsTable.M==1);
if sum(pos)>0
    plot(segCoord(pos,:)',log2(ones(2,1)*2'+1),'k','linewidth',4);
end
%plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)','k','linewidth',4);
ticks=[0 2.^[0:7]];
tickpos=log2(ticks+1);
set(gca,'YTick',tickpos,'YTickLabel',ticks,'tickDir','out');
set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrList,'FontSize',6);
axis([0 max(segCoord(:,2)) log2(1) log2(max(segsTable.N)+1)])
title('Copy Number','FontSize',10);

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 7.5 10]);
print(gcf,'-dpng',[inputParam.outName '.cnaPlot.png'],'-r300');
close(gcf);


%%% plot het AF
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1])
    T=Tcell{j};
    bIdx=T.ApopAF>=T.BpopAF;
    aIdx=T.ApopAF<T.BpopAF;
    AF(bIdx)=(T.BCountF(bIdx)+T.BCountR(bIdx))./T.ReadDepthPass(bIdx);
    AF(aIdx)=(T.ACountF(aIdx)+T.ACountR(aIdx))./T.ReadDepthPass(aIdx);
    for i=1:2:maxChr
        scatter(Tcoord(hetPos & T.Chr==i),AF(hetPos & T.Chr==i),1,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
        hold on;
    end
    for i=2:2:maxChr
        scatter(Tcoord(hetPos & T.Chr==i),AF(hetPos & T.Chr==i),1,'.','MarkerFaceColor',[0.8235    0.7059    0.5490],'MarkerEdgeColor',[0.8235    0.7059    0.5490])
    end
end

%%% plot somatic AF
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1])
    hold on;
    T=Tcell{j};
    bIdx=T.ApopAF>=T.BpopAF;
    aIdx=T.ApopAF<T.BpopAF;
    AF(bIdx)=(T.BCountF(bIdx)+T.BCountR(bIdx))./T.ReadDepthPass(bIdx);
    AF(aIdx)=(T.ACountF(aIdx)+T.ACountR(aIdx))./T.ReadDepthPass(aIdx);
    for i=1:size(f,2)
        pos=cloneId(:,1)==i & somPos;
        if(rem(i,2)==0)
            scatter(Tcoord(pos),AF(pos),5,'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',[0 0 0]);
        else
            scatter(Tcoord(pos),AF(pos),5,'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
        end
    end
    axis([0 max(segCoord(:,2)) 0 1]);
    set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrList,'tickDir','out','FontSize',8);
    title('Allele Frequencies','FontSize',10);
end

for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,j*5);
    hold on;
    for i=1:size(f,2)
        b=bar(i,f(j,i));
        b.FaceColor=colors(i,:);
        if(rem(i,2)==0)
            b.EdgeColor=[0 0 0];
        else
            b.EdgeColor=[1 1 1];
        end
    end
    ylim([0 1]);
    xlim([0 size(f,2)+1]);
    ylabel('Sample Fraction','FontSize',10);
end

for j=1:length(Tcell)
    cnCorr(:,j)=segsTable.F(:,j).*segsTable.M./segsTable.N+(1-segsTable.F(:,j))*0.5;
    cnCorr(segsTable.N==0,j)=0.5;
end

for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1])
    for i=1:size(f,2)
        pos=segsTable.cnaIdx==i & (segsTable.N~=2 | segsTable.M~=1);
        plot(segCoord(pos,:)',ones(2,1)*cnCorr(pos,j)','Color','k','linewidth',2);
        plot(segCoord(pos,:)',1-ones(2,1)*cnCorr(pos,j)','Color','k','linewidth',2);
        %plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)',colors(i),'linewidth',4);
    end
    pos=(segsTable.N==2 & segsTable.M==1);
    plot(segCoord(pos,:)',ones(2,1)*cnCorr(pos,j)','k','linewidth',2);
    %plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)','k','linewidth',4);
    %ticks=[0 2.^[1:7]];
    %tickpos=log2(median(segsTable(:,7)).*ticks/2+(1-median(segsTable(:,7))));
    %set(gca,'YTick',tickpos,'YTickLabel',ticks,'tickDir','out');
    set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrList,'FontSize',6);
    axis([0 max(segCoord(:,2)) 0 1])
    title(samples{j},'FontSize',10,'Interpreter','none');
end
subplot(length(Tcell)+1,5,[5*(length(Tcell)+1)-4 5*(length(Tcell)+1)-1])
hold on;
for i=1:size(f,2)
    pos=segsTable.cnaIdx==i & (segsTable.N~=2 | segsTable.M~=1);
    if(rem(i,2)==0)
        plot(segCoord(pos,:)',ones(2,1)*(segsTable.M(pos)'./segsTable.N(pos)'),'Color',[0 0 0],'linewidth',3);
        plot(segCoord(pos,:)',ones(2,1)*(segsTable.M(pos)'./segsTable.N(pos)'),'Color',colors(i,:),'linewidth',2);
    else
        plot(segCoord(pos,:)',ones(2,1)*(segsTable.M(pos)'./segsTable.N(pos)'),'Color',colors(i,:),'linewidth',3);
    end
    %plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)',colors(i),'linewidth',4);
end
pos=(segsTable.N==2 & segsTable.M==1);
plot(segCoord(pos,:)',ones(2,1)*(segsTable.M(pos)'./segsTable.N(pos)'),'k','linewidth',2);
%plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)','k','linewidth',4);
%ticks=[0 2.^[0:7]];
%tickpos=log2(ticks+1);
%set(gca,'YTick',tickpos,'YTickLabel',ticks,'tickDir','out');
set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrList,'FontSize',6);
axis([0 max(segCoord(:,2)) 0 0.5])
title('M/N','FontSize',10);

%%% print plot
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 10 7.5]);
print(gcf,'-dpng',[inputParam.outName '.vafPlot.png'],'-r300');
close(gcf);
return;