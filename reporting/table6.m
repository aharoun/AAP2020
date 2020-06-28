%% table 6

tb6.entry         = (eqINCF.z/eqIN.z-1)*100;
tb6.x             = (eqINCF.xStar(1:5)'./eqIN.xStar(1:5)' - 1)*100;
tb6.xAverage      = ((eqINCF.muHigh'*eqINCF.xStar)/(eqIN.muHigh'*eqIN.xStar) - 1)*100;
tb6.aveFirmSize   = ((eqIN.hM + eqIN.lM)./(eqINCF.hM + eqINCF.lM) - 1)*100;
tb6.shareHighFirm = ((eqINCF.hM/(eqINCF.hM + eqINCF.lM))./(eqIN.hM/(eqIN.hM + eqIN.lM)) - 1)*100;
tb6.empShareSmall = ((eqINCF.fAll(1)*eqINCF.emp(1)/eqINCF.Ltotal)/(eqIN.fAll(1)*eqIN.emp(1)/eqIN.Ltotal) - 1)*100;
tb6.tauChange     = (eqINCF.tau/eqIN.tau - 1)*100;
tb6.meanEmp       = (mINCF.meanEmp     ./mIN.meanEmp - 1)*100;
tb6.shareSmall    = (mINCF.shareOneProd./mIN.shareOneProd - 1)*100;
tb6.managerShare  = (mINCF.shareMan./mIN.shareMan - 1)*100;


fd = fopen(['output' filesep 'Table6.txt'],'w+');

disp('Table 6 : Increasing the Delegation Efficiency in India: Firm-level implications');
fprintf(fd,'%30s\n\n','Table 6:Increasing the Delegation Efficiency in India: Firm-level implications');
fprintf(fd,'------------------------------------------------------------------------------------------------------\n');

fprintf(fd,'%27s\t %10s\t %10s\t %10s\t %10s\t %10s\t %10s\n',' ','Average','n=1','n=2','n=3','n=4','n=5');
fprintf(fd,'%27s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','Expansion rate',tb6.xAverage,tb6.x);

fprintf(fd,'------------------------------------------------------------------------------------------------------\n');
fprintf(fd, '%27s\t %10.2f\n','Entry intensity',tb6.entry);
fprintf(fd, '%27s\t %10.2f\n','Creative destruction',tb6.tauChange);
fprintf(fd, '%27s\t %10.2f\n','Share of manager',tb6.managerShare);
fprintf(fd, '%27s\t %10.2f\n','Average firm size',tb6.aveFirmSize);
fprintf(fd, '%27s\t %10.2f\n','Share of high type firms',tb6.shareHighFirm);
fprintf(fd, '%27s\t %10.2f\n','Empl. share of small firms',tb6.empShareSmall);
fprintf(fd,'------------------------------------------------------------------------------------------------------\n');

fprintf(fd,'%27s\t %10s\t %10s\t %10s\t %10s\t %10s\t %10s\n',' ','5<=','6-10','11-15','16-20','21-25','+26');
fprintf(fd,'%27s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','Average firm size',tb6.meanEmp);
fprintf(fd,'%27s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','Share of small firms',tb6.shareSmall);

fprintf(fd,'------------------------------------------------------------------------------------------------------\n');
fclose(fd);
