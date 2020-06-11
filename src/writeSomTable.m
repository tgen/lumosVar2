function writeSomTable(Tcell,somPos,cloneId,sampleFrac,somaticDetected,outName)

somTable=Tcell{1}(somPos,{'Chr','Pos'});
somTable.Ref=int2ntIndels(Tcell{1}.RefComb(somPos));
somTable.A_allele=int2ntIndels(Tcell{1}.Acomb(somPos));
somTable.B_allele=int2ntIndels(Tcell{1}.Bcomb(somPos));
somTable.cloneGroup=cloneId(somPos,1);
somTable.detectStr=num2str(somaticDetected(somPos,:),'%-d');
somTable.sampleFrac=sampleFrac(somPos,:);
for i=1:length(Tcell)
	somTable.Depth(:,i)=Tcell{i}.ReadDepthPass(somPos);
    somTable.Bcounts(:,i)=Tcell{i}.BcountsComb(somPos);
end
somTable.N=Tcell{1}.NumCopies(somPos);
somTable.M=Tcell{1}.MinAlCopies(somPos);


writetable(somTable,[outName '.somaticPass.txt'],'Delimiter','\t','FileType','text','WriteVariableNames',1);
