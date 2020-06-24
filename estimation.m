%% Estimation routine

clear all
close all
clc

global scoreBest

scoreBest = Inf;
vers = '';  

if any(strfind(vers,'Payment')) 
	parinit = ones(1,14);  % sigma is estimated for countries separately
else
	parinit = ones(1,13);  % sigma is common across countries
end

mopts = optimset('Display','iter','MaxFunEvals',1500,'MaxIter',1500);
[parfin,~] = fminsearch(@(x)estimationObj(x,vers),parinit,mopts);
% final result
[score,pIN,eqIN,mIN,pUS,eqUS,mUS,meanReg] = estimationObj(parfin,vers);

