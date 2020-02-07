function plotCNAandVAF(exonRD,segsTable,CNAscale,fIn,Tcell,somPos,hetPos,cloneId,inputParam)
%plotCNAandVAF - make copy number and VAF plots
%
% Syntax:  plotCNAandVAF(exonRD,segsTable,CNAscale,fIn,Tcell,somPos,hetPos,cloneId,inputParam)
%
% Inputs:
%   exonRD: cell array of matrices of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segsTable: table of segment data with variables: {'Chr','StartPos','EndPos',
%       'N','M','F','cnaIdx'}
%   CNAscale: copy number centering parameter
%   fIn: sample fraction matrix
%   Tcell: cell array of tables with columns: {'Chr','Pos','ReadDepthPass'}
%   somPos: logical vector indicating positions called somatic
%   hetPos: logical vector of positions called heterozygous
%   cloneId: index of clonal variant group
%   inputParam: structure of parameters
%
% Outputs:
%    output file names/paths specified by inputParam.outName
%    *.cnaPlot.pdf - copy number plot
%    *.vafPlot.pdf - plots of variant allele fractions
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain, callSNVmulti, callVariants

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 9-May-2018
%------------- BEGIN CODE --------------

%%%fill in sample fractions for normal sample
if inputParam.NormalSample>0
    f=[zeros(length(Tcell),inputParam.numClones) ones(length(Tcell),1)];
    tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
    f(tIdx,1:end-1)=[fIn];
else
   f=[fIn ones(length(Tcell),1)];
end

%%% transform chromosome coord to linear coord
chrTable=inputParam.chrTable;
chrNum=chrTable.chrIdx;
T=Tcell{1};
for i=1:length(chrNum)
    chrLen(i)=max(segsTable{segsTable.Chr==chrNum(i),3})+1E7;
end
chrOffset=[0; cumsum(chrLen(1:end-1))'];
for i=1:length(chrNum)
    segCoord(segsTable.Chr==chrNum(i),1)=(segsTable{segsTable.Chr==chrNum(i),2}+chrOffset(i))/1E6;
    segCoord(segsTable.Chr==chrNum(i),2)=(segsTable{segsTable.Chr==chrNum(i),3}+chrOffset(i))/1E6;
    exonCoord(exonRD{1}(:,1)==chrNum(i),1)=(exonRD{1}(exonRD{1}(:,1)==chrNum(i),2)+chrOffset(i))/1E6;
    exonCoord(exonRD{1}(:,1)==chrNum(i),2)=(exonRD{1}(exonRD{1}(:,1)==chrNum(i),3)+chrOffset(i))/1E6;
    Tcoord(T.Chr==chrNum(i),1)=(T.Pos(T.Chr==chrNum(i))+chrOffset(i))/1E6;
end

%%%get chromosomes


%%% calculate log2FC
for i=1:length(Tcell)
    log2FC(:,i)=log2((CNAscale(i)./2).*exonRD{i}(:,4)./exonRD{i}(:,5));
    Nlog2R(:,i)=log2(segsTable.F(:,i).*segsTable.N/2+(1-segsTable.F(:,i)));
end

%%% plot exon log2FC
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1]);
    for k=1:2:length(chrNum)
	i=chrNum(k);
        idx=exonRD{j}(:,1)==i;
        scatter(mean(exonCoord(idx,:),2),log2FC(idx,j),1,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
        hold on;
    end
    for k=2:2:length(chrNum)
        i=chrNum(k);
        idx=exonRD{j}(:,1)==i;
        scatter(mean(exonCoord(idx,:),2),log2FC(idx,j),1,'.','MarkerFaceColor',[0.8235    0.7059    0.5490],'MarkerEdgeColor',[0.8235    0.7059    0.5490])
        hold on;
    end
    ylabel('log2(FoldChange)','FontSize',10);
end

%%%make sample fraction bar plots
colors=linspecer(size(f,2)-1);
colors=[colors; [0 0 0]];
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,j*5);
    hold on;
    for i=1:size(f,2)-1
        b=bar(i,f(j,i));
        b.FaceColor=colors(i,:);
    end
    axis tight;
    ylim([0 1]);
    ylabel('Sample Fraction','FontSize',10);
end

%%% plot fitted segments
samples=regexp(inputParam.sampleNames,',','split');
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1])
    plot(reshape(segCoord',[],1),reshape(ones(2,1)*Nlog2R(:,j)',[],1),'k','linewidth',1);
    hold on;
    for i=1:size(f,2)
        pos=segsTable.cnaIdx==i & (segsTable.N~=2 | segsTable.M~=1);
        plot(segCoord(pos,:)',ones(2,1)*Nlog2R(pos,j)','Color',colors(i,:),'linewidth',4); 
    end
    pos=(segsTable.N==2 & segsTable.M==1);
    plot(segCoord(pos,:)',ones(2,1)*Nlog2R(pos,j)','k','linewidth',3);
    set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrTable.chrName,'FontSize',6);
    axis([0 max(segCoord(:,2)) min(Nlog2R(:))-1 max(Nlog2R(:))+1])
    title(samples{j},'FontSize',10,'Interpreter','none');
end

%%%plot integer copy number
subplot(length(Tcell)+1,5,[5*(length(Tcell)+1)-4 5*(length(Tcell)+1)-1])
hold on;
plot(reshape(segCoord',[],1),reshape(log2(ones(2,1)*segsTable.N'+1),[],1),'k','linewidth',1);
for i=1:size(f,2)
    pos=segsTable.cnaIdx(:,1)==i & (segsTable.N~=2 | segsTable.M~=1);
    plot(segCoord(pos,:)',log2(ones(2,1)*segsTable.N(pos)'+1),'Color',colors(i,:),'linewidth',3);
end
pos=(segsTable.N==2 & segsTable.M==1);
if sum(pos)>0
    plot(segCoord(pos,:)',log2(ones(2,1)*2'+1),'k','linewidth',4);
end
ticks=[0 2.^[0:7]];
tickpos=log2(ticks+1);
set(gca,'YTick',tickpos,'YTickLabel',ticks,'tickDir','out');
set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrTable.chrName,'FontSize',6);
axis([0 max(segCoord(:,2)) log2(1) log2(max(segsTable.N)+1)])
title('Copy Number','FontSize',10);

%%%print copy number plot
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 7.5 10]);
print(gcf,'-dpdf',[inputParam.outName '.cnaPlot.pdf'],'-r300');
close(gcf);


%%% plot het AF
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1])
    T=Tcell{j};
    bIdx=T.ApopAFcomb>=T.BpopAFcomb;
    aIdx=T.ApopAFcomb<T.BpopAFcomb;
    AF(bIdx)=(T.BCountF(bIdx)+T.BCountR(bIdx))./T.ReadDepthPass(bIdx);
    AF(aIdx)=(T.ACountF(aIdx)+T.ACountR(aIdx))./T.ReadDepthPass(aIdx);
    for k=1:2:length(chrNum)
        i=chrNum(k);
        scatter(Tcoord(hetPos & T.Chr==i),AF(hetPos & T.Chr==i),1,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
        hold on;
    end
    for k=2:2:length(chrNum)
        i=chrNum(k);
        scatter(Tcoord(hetPos & T.Chr==i),AF(hetPos & T.Chr==i),1,'.','MarkerFaceColor',[0.8235    0.7059    0.5490],'MarkerEdgeColor',[0.8235    0.7059    0.5490])
    end
end

%%% plot somatic AF
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1])
    hold on;
    T=Tcell{j};
    bIdx=T.ApopAFcomb>=T.BpopAFcomb;
    aIdx=T.ApopAFcomb<T.BpopAFcomb;
    AF(bIdx)=(T.BcountsComb(bIdx))./T.ReadDepthPass(bIdx);
    AF(aIdx)=(T.AcountsComb(aIdx))./T.ReadDepthPass(aIdx);
    for i=1:size(f,2)
        pos=cloneId(:,1)==i & somPos;
        scatter(Tcoord(pos),AF(pos),5,'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
    end
    axis([0 max(segCoord(:,2)) 0 1]);
    set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrTable.chrName,'tickDir','out','FontSize',8);
    title('Allele Frequencies','FontSize',10);
end

%%%make sample fraction bar plots
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,j*5);
    hold on;
    for i=1:size(f,2)-1
        b=bar(i,f(j,i));
        b.FaceColor=colors(i,:);
    end
    axis tight;
    ylim([0 1]);
    ylabel('Sample Fraction','FontSize',10);
end

%%%find expected hetAF per segment
for j=1:length(Tcell)
    cnCorr(:,j)=(segsTable.F(:,j).*segsTable.M+(1-segsTable.F(:,j)))./(segsTable.F(:,j).*segsTable.N+2*(1-segsTable.F(:,j)));
    cnCorr(segsTable.N==0,j)=0.5;
end

%%% plot segment expected het AF
for j=1:length(Tcell)
    subplot(length(Tcell)+1,5,[j*5-4 j*5-1])
    for i=1:size(f,2)
        pos=segsTable.cnaIdx==i & (segsTable.N~=2 | segsTable.M~=1);
        plot(segCoord(pos,:)',ones(2,1)*cnCorr(pos,j)','Color','k','linewidth',2);
        plot(segCoord(pos,:)',1-ones(2,1)*cnCorr(pos,j)','Color','k','linewidth',2);
     end
    pos=(segsTable.N==2 & segsTable.M==1);
    plot(segCoord(pos,:)',ones(2,1)*cnCorr(pos,j)','k','linewidth',2);
    set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrTable.chrName,'FontSize',6);
    axis([0 max(segCoord(:,2)) 0 1])
    title(samples{j},'FontSize',10,'Interpreter','none');
end

%%% plot M/N (expected hetAF for 100% tumor)
subplot(length(Tcell)+1,5,[5*(length(Tcell)+1)-4 5*(length(Tcell)+1)-1])
hold on;
for i=1:size(f,2)
    pos=segsTable.cnaIdx==i & (segsTable.N~=2 | segsTable.M~=1);
    plot(segCoord(pos,:)',ones(2,1)*(segsTable.M(pos)'./segsTable.N(pos)'),'Color',colors(i,:),'linewidth',3);
end
pos=(segsTable.N==2 & segsTable.M==1);
plot(segCoord(pos,:)',ones(2,1)*(segsTable.M(pos)'./segsTable.N(pos)'),'k','linewidth',2);
set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',chrTable.chrName,'FontSize',6);
axis([0 max(segCoord(:,2)) 0 0.5])
title('M/N','FontSize',10);

%%% print allele fraction plot
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 7.5 10]);
print(gcf,'-dpdf',[inputParam.outName '.vafPlot.pdf'],'-r300');
close(gcf);
return;
