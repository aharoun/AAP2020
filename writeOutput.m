function [] = writeOutput(param,m,filename,title)
global alg 
if isstruct(param)
   p = cell2mat(struct2cell(param));    
else
   p = param;
end

[~,paramNames]    = readParam(alg.paramsFile);
% Printing the final results
% construct cell array
result      = cell(m.nMoment,6);
result(:,1) = num2cell(m.momentModel);
result(:,2) = num2cell(m.momentData);
result(:,3) = num2cell(1:m.nMoment);
result(1:length(m.momentName),4) = m.momentName;
result(:,5) = num2cell(m.momentWgt);
result(:,6) = num2cell(m.mErr);

fd = fopen([filename '.txt'],'w+');
fprintf(fd,'%10s\n\n\n',title);
fprintf(fd,'%10s\n','MOMENTS');
fprintf(fd,'%10s\t %10s\t %5s\t %30s\t %10s\n','model','data','#','description','weight');

for i=1:length(m.momentName)
   fprintf(fd,'%10.6f\t %10.6f\t %5i\t %30s\t %10.6f\t\n',result{i,1},result{i,2},result{i,3},result{i,4},result{i,5});
end
fprintf(fd,'\n\n');
fprintf(fd,'%10s\n','PARAMETERS');
  
for i = 1: length(p)
   fprintf(fd,'%8s:\t ',char(paramNames{i}));
   fprintf(fd,'%2.15f\n',p(i));
end
fprintf(fd,'\n\n');
fprintf(fd,'Country:\t %10s\n',alg.cnt);
fprintf(fd,'\n\n');
fclose(fd);
  
end
