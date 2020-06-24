%% table 5

fd = fopen(['output' filesep 'Table5.txt'],'w+');

disp('Table 5 : Creative Destruction and Selection in India and the U.S.');
fprintf(fd,'%50s\n\n\n','Table 5:  Creative Destruction and Selection in India and the U.S.');
fprintf(fd,'%50s\t %10s\t %10s\n',' ','India', 'US');
fprintf(fd,'------------------------------------------------------------------------------------\n');

fprintf(fd,'%50s\t %10.4f\t %10.4f\t\n','Rate of creative destruction',eqIN.tau,eqUS.tau);
fprintf(fd,'%50s\t %10.4f\t %10.4f\t\n','Share of high-type firms upon entry',pIN.delta,pUS.delta);
fprintf(fd,'%50s\t %10.4f\t %10.4f\t\n','Long-run share of high-type firms',eqIN.hM/(eqIN.hM + eqIN.lM),eqUS.hM/(eqUS.hM + eqUS.lM));
fprintf(fd,'%50s\t %10.4f\t %10.4f\t\n','Long-run employment share of high-type firms',eqIN.fHigh'*eqIN.emp,eqUS.fHigh'*eqUS.emp);
fprintf(fd,'%50s\t %10.4f\t %10.4f\t\n','Long-run share of high-type firms among age 21-25',mIN.shareHighAgeBin1(5),mUS.shareHighAgeBin1(5));

fprintf(fd,'------------------------------------------------------------------------------------\n');
fclose(fd);
