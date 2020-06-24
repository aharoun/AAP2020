%% table 9
robustnessList{1} = 'Baseline';

fd = fopen(['output' filesep 'Table9.txt'],'w+');
disp('Table 9 : Robustness')
fprintf(fd,'%50s\n\n','Table 9: Robustness');
fprintf(fd,'------------------------------------------------------------------------------------------------------------\n');
fprintf(fd,'%16s\t %10s\t %10s\t %10s\t %10s\t %10s\t %10s\t %10s\n','','tauIN','tauUS','dTauIN','dQIN2US','dYIN2US','dAveFirmSize','dShareSmall2125');
for i = 1:length(robustnessList)
	fprintf(fd,'%16s\t %10.3f\t %10.3f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n',...
	robustnessList{i},rb.tauIN(i),rb.tauUS(i),rb.dTauIN(i),rb.dQIN2US(i),rb.dYIN2US(i),rb.dAverageFirmSize(i),rb.dShareOneProd2125(i));
end

fprintf(fd,'------------------------------------------------------------------------------------------------------------\n');
fclose(fd);
