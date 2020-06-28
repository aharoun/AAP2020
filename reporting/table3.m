% Table 3: Parameters for the U.S. and India

fd = fopen(['output' filesep 'Table3.txt'],'w+');
disp('Table 3 : Estimated Parameters for the U.S. and India')
fprintf(fd,'%10s\n\n','Table 3: Estimated Parameters for the U.S. and India');
fprintf(fd,'%10s\t %10s\t %10s\n','Parameter','US','India');
fprintf(fd,'---------------------------------------\n');

paramNames = fieldnames(pUS);
for i=1:length(paramNames)-1
   fprintf(fd,'%10s\t %10.5f\t %10.5f\t\n', paramNames{i},pUS.(paramNames{i}),pIN.(paramNames{i}));
end

fprintf(fd,'---------------------------------------\n');
fclose(fd);

