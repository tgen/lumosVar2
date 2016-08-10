function writeSomTable(Tcell,filtPos,P)

posList=[];
for j=1:size(Tcell,2)
    somPos=P{i,j}.Somatic>inputParam.pSomaticThresh & filtPos{j} & min([T.ApopAF T.BpopAF],[],2)<inputParam.maxSomPopFreq;
    posList=[posList; T{somPos,1:2}];
end
posList=unique(posList,'rows');
somTable=array2table(posList,'VariableNames',{'Chr','Pos'});

for i=1:size(Tcell,2)
    T=Tcell{i};
    aIdx=T.ApopAF>=T.BpopAF;
    AF(aIdx)=(T.BCountF(aIdx)+T.BCountR(aIdx))./T.ReadDepthPass(aIdx);
    bIdx=T.ApopAF<T.BpopAF;
    AF(bIdx)=(T.ACountF(bIdx)+T.ACountR(bIdx))./T.ReadDepthPass(bIdx);
    sampleFrac=AF.*(T.cnaF.*T.NumCopies+2*(1-T.cnaF))./(T.cnaF.*(T.NumCopies-T.MinAlCopies)+(1-T.cnaF));
    %sampleFrac=AF.*(T.cnaF.*T.NumCopies+2*(1-T.cnaF))./(T.NumCopies-T.MinAlCopies);
    [lia,locb]=ismember(posList,[T.Chr T.Pos],'Rows');
    somTable.pTrust(lia,i)=P{1,i}.postTrust(locb(lia),2);
    somTable.pArtifact(lia,i)=P{1,i}.postArtifact(locb(lia),2);
    somTable.pSomatic(lia,i)=P{end,i}.Somatic(locb(lia));
    somTable.AF(lia,i)=AF(locb(lia));
    somTable.clonalFrac(lia,i)=f{end,i}(cloneId{end,i}(locb))';
    somTable.sampleFrac(lia,i)=sampleFrac(locb(lia));
end
    
    
somTable.sampleFrac(somTable.pSomatic<0.8)=0;   

filtIdx=find(min(somTable.pArtifact,[],2)>0.99);
sampleFracFilt=somTable.sampleFrac(min(somTable.pArtifact,[],2)>0.99,:);
clustId=clusterdata(sampleFracFilt,'maxclust',9,'linkage','average','distance','cityblock');

clonalFracFilt=somTable.clonalFrac(min(somTable.pArtifact,[],2)>0.99,:);
clustId=clusterdata(clonalFracFilt,'maxclust',16,'linkage','average');

for i=1:16
    subplot(4,4,i)
    plot(sampleFracFilt(i==clustId,:)');
end
for i=1:4
    sampleFracFilt=sortrows(sampleFracFilt,5-i);
end

pcolor([sampleFracFilt ones(size(sampleFracFilt,1),1); ones(1,size(sampleFracFilt,2)+1)]);