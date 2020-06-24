%% table 7

tb7.entry         = (eqUSCF.z/eqUS.z-1)*100;
tb7.xAverage      = ((eqUSCF.muHigh'*eqUSCF.xStar)/(eqUS.muHigh'*eqUS.xStar) - 1)*100;
tb7.aveFirmSize   = ((eqUS.hM + eqUS.lM)./(eqUSCF.hM + eqUSCF.lM) - 1)*100;
tb7.shareHighFirm = ((eqUSCF.hM/(eqUSCF.hM + eqUSCF.lM))./(eqUS.hM/(eqUS.hM + eqUS.lM)) - 1)*100;
tb7.empShareSmall = ((eqUSCF.fAll(1)*eqUSCF.emp(1)/eqUSCF.Ltotal)/(eqUS.fAll(1)*eqUS.emp(1)/eqUS.Ltotal) - 1)*100;
tb7.tauChange     = (eqUSCF.tau/eqUS.tau - 1)*100;
tb7.meanEmp       = (mUSCF.meanEmp     ./mUS.meanEmp - 1)*100;
tb7.shareSmall    = (mUSCF.shareOneProd./mUS.shareOneProd - 1)*100;
tb7.managerShare  = (mUSCF.shareMan./mUS.shareMan - 1)*100;


fd = fopen(['output' filesep 'Table7.txt'],'w+');

disp('Table 7 : Decreasing Delegation Efficiency in the U.S.');
fprintf(fd,'%30s\n\n','Table 7: Decreasing Delegation Efficiency in the U.S.');
fprintf(fd,'------------------------------------------\n');
fprintf(fd, '%27s\t %10.2f\n','Average expansion',tb7.xAverage);
fprintf(fd, '%27s\t %10.2f\n','Entry intensity',tb7.entry);
fprintf(fd, '%27s\t %10.2f\n','Creative destruction',tb7.tauChange);
fprintf(fd, '%27s\t %10.2f\n','Average firm size',tb7.aveFirmSize);
fprintf(fd, '%27s\t %10.2f\n','Share of high type firms',tb7.shareHighFirm);
fprintf(fd, '%27s\t %10.2f\n','Empl. share of small firms',tb7.empShareSmall);
fprintf(fd, '%27s\t %10.2f\n','Share of manager',tb7.managerShare);
fprintf(fd,'------------------------------------------\n');
fclose(fd);
