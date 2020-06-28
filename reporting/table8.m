%% Table 8
% Construct output path
yUSInitBase = 1.0; 
qUSt0       = yUSInitBase/(eqTransUSBase.lProdS(1)*eqTransUSBase.Mt(1));
qINt0       = qUSt0/qUS2qIndt0;
yINInitBase = eqIN.lProdS*eqIN.M*qINt0;
yINCFInit   = eqTransINCF. lProdS(1)*eqTransINCF. Mt(1)*qINt0;

yPathUSBase = cumprod([yUSInitBase; (1 + eqTransUSBase. gY(1:end-1)*alg.dt)]);
yPathINBase = cumprod([yINInitBase; (1 + eqTransINBase. gY(1:end-1)*alg.dt)]);
yPathINCF   = cumprod([yINCFInit;(1 + eqTransINCF.gY(1:end-1)*alg.dt)]); 

% Construcy Q path
qPathUSBase = cumprod([qUSt0; (1 + eqTransUSBase. gQ(1:end-1)*alg.dt)]);
qPathINBase = cumprod([qINt0; (1 + eqTransINBase. gQ(1:end-1)*alg.dt)]);
qPathINCF   = cumprod([qINt0;(1 + eqTransINCF.gQ(1:end-1)*alg.dt)]); 

% Construct consumption path
cPathUSBase = yPathUSBase.* (1 - eqTransUSBase. RDCost - eqTransUSBase. entryCost);
cPathINBase = yPathINBase.* (1 - eqTransINBase. RDCost - eqTransINBase. entryCost);
cPathINCF   = yPathINCF.*(1 - eqTransINCF.RDCost - eqTransINCF.entryCost);
%--------------------------------------------------------------------------------------------------------------------------------

gIN2000 = log(pIN.gamma*real(qUS2qINPathModel(timePath+1985 ==2000)).^pIN.lambda)*eqIN.tau;

qRatio2000 = 1/qUS2qINPathModel(timePath+1985 ==2000);
yRatio2000 = (eqIN.lProdS*eqIN.M)/(eqUS.lProdS*eqUS.M)*qRatio2000;
cRatio2000 = yRatio2000*(1.0 - eqTransINBase.RDCost(1) - eqTransINBase.entryCost(1))/(1.0 - eqTransUSBase.RDCost(1) - eqTransUSBase.entryCost(1));

t = alg.timePath + 2010-alg.dt;

y = [2010, 2020, 2030, t(end)];

result = [gIN2000 gIN2000 qRatio2000 qRatio2000 yRatio2000 yRatio2000 cRatio2000 cRatio2000];
for i = 1:length(y)
	result = [result; eqTransINBase.gQ(t==y(i)) eqTransINCF.gQ(t==y(i))...
					  qPathINBase(t==y(i))/qPathUSBase(t==y(i))...
	                  qPathINCF(t==y(i))/qPathUSBase(t==y(i))... 
	                  yPathINBase(t==y(i))/yPathUSBase(t==y(i))...
	                  yPathINCF(t==y(i))/yPathUSBase(t==y(i))...
	                  cPathINBase(t==y(i))/cPathUSBase(t==y(i))...
	                  cPathINCF(t==y(i))/cPathUSBase(t==y(i))]; 
end

result = result'*100;


fd = fopen(['output' filesep 'Table8.txt'],'w+');

disp('Table 8 : Increasing Delegation Efficiency in India: Macroeconomic Implications');
fprintf(fd,'%50s\n\n','Table 8: Increasing Delegation Efficiency in India: Macroeconomic Implications');
fprintf(fd,'------------------------------------------------------------------------------------\n');
fprintf(fd,'%10s\t %10s\t %10s\t %10s\t %10s\t %10s\n',' ','2000','2010','2020','2030','S.S.');
fprintf(fd,'------------------------------------------------------------------------------------\n');
fprintf(fd,'%60s\n','Productivity Growth');
fprintf(fd,'%10s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','Baseline',result(1,:));
fprintf(fd,'%10s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','CF',result(2,:));
fprintf(fd,'------------------------------------------------------------------------------------\n');
fprintf(fd,'%60s\n','Relative productivity ');
fprintf(fd,'%10s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','Baseline',result(3,:));
fprintf(fd,'%10s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','CF',result(4,:));
fprintf(fd,'------------------------------------------------------------------------------------\n');
fprintf(fd,'%60s\n','Relative income pc');
fprintf(fd,'%10s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','Baseline',result(5,:));
fprintf(fd,'%10s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','CF',result(6,:));
fprintf(fd,'------------------------------------------------------------------------------------\n');
fprintf(fd,'%60s\n','Relative Consumption ');
fprintf(fd,'%10s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','Baseline',result(7,:));
fprintf(fd,'%10s\t %10.2f\t %10.2f\t %10.2f\t %10.2f\t %10.2f\n','CF',result(8,:));

fprintf(fd,'------------------------------------------------------------------------------------\n');
fclose(fd);
