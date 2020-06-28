% Table 4: Moments for the U.S. and India

fd = fopen(['output' filesep 'Table4.txt'],'w+');

disp('Table 4 : Estimated Parameters for the U.S. and India');
fprintf(fd,'%19s\n\n\n','Table 4:  Estimated Parameters for the U.S. and India');
fprintf(fd,'%19s\t %10s\t %10s\t %10s\t %10s\n','Moment','US-data','US-model','India-data','India-model');
fprintf(fd,'--------------------------------------------------------------------\n');

for i=1:length(mUS.momentName)
   fprintf(fd,'%19s\t %10.4f\t %10.4f\t %10.4f\t %10.4f\t\n',mUS.momentName{i}, mUS.momentData(i),mUS.momentModel(i), mIN.momentData(i),mIN.momentModel(i));
end

fprintf(fd,'%19s\t %10s\t %10s\t %10.4f\t %10.4f\t\n','Treatment eff.', 'N/A','N/A', mIN.data.bloomMoment,bloomMomentModel);
fprintf(fd,'%19s\t %10s\t %10s\t %10.4f\t %10.4f\t\n','relManShareMigrant', 'N/A','N/A', mIN.data.shareManMigrantUS/mIN.data.shareManMigrantIN,...
                                                                            (pIN.muM/pUS.muM)^pIN.vartheta*(mUS.data.shareMan/mIN.data.shareMan) );
fprintf(fd,'%19s\t %10.4f\t %10.4f\t %10.4f\t %10.4f\t\n','varLogManWage',  mUS.data.varLogomegaM,mUS.varLogomegaM , mIN.data.varLogomegaM,mIN.varLogomegaM);

fprintf(fd,'--------------------------------------------------------------------\n');
fclose(fd);
