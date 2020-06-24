function [] = writeParam(values,names,fname)
% Writes parameters to file

if nargin<3
    error('Not enough input arguments')
end

if isstruct(values)
    values = cell2mat(struct2cell(values));
end


if length(values)~=length(names)
    error('Length of names and values does not match')
end


fid = fopen(fname,'w+');
  
for i = 1: length(values)-1
    fprintf(fid,'%8s : ',char(names{i}));
    fprintf(fid,'%2.15f\n',values(i));
end

fprintf(fid,'%8s : ',char(names{i+1}));
fprintf(fid,'%2.15f',values(i+1));
  
fclose(fid);
  
end