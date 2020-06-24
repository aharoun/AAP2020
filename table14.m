%% table 14

disp('Table 14: Moment Sensitivity')
fd = fopen(['output' filesep 'Table14.txt'],'w+');
fprintf(fd,'%50s\n\n','Table 14: Moment Sensitivity');
fprintf(fd,['%14s\t ' repmat('%8s\t ',1,7), '\n'],momSensitivityMatrix{1,:});
fprintf(fd,'-------------------------------------------------------------------------------------------------\n');
for i = 2:8
	fprintf(fd,['%14s\t ' repmat('%8.3f\t ',1,7), '\n'],momSensitivityMatrix{i,:});
end

fprintf(fd,'--------------------------------------------------------------------------------------------------');
fclose(fd);
