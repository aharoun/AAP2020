%% REPLICATION FILE

% Replication file for "Lack of Selection and Limits to Delegation: Firm Dynamics in Developing Countries" 
% by Ufuk Akcigit, Harun Alp and Michael Peters
% April 2020

clear all;
close all;
clc;
global alg

% -------------------------------------------------------------------------
% ------------------------- BASELINE --------------------------------------
% -------------------------------------------------------------------------
disp('------------------');
disp(' BASELINE RESULTS ')
disp('------------------');
[eqUS,mUS,pUS] = solver('US');					  % Solve and get moments
[eqIN,mIN,pIN] = solver('India');				  % Solve and get moments
bloomMomentModel = bloomExercise(eqIN,pIN,eqUS);  % Bloom regression moment
%---------------------------------------------------------------------------

disp('Solving India Transition')
[~,qUS2qINPathModel,timePath] = calibrateDiffusion(pIN,eqIN,pUS,eqUS);
qUS2qIndt0 = qUS2qINPathModel(timePath+1985 ==2010);
pIN.qUS2qIndInit = qUS2qIndt0;  % initial aggregate Q ratios

initAlg('India');
[eqSS0,p0] = solveBGP(pIN,'India'); % Initial steady state
[eqSS1,p1] = solveBGP(pIN,'India'); % Final steady state

[initState,initSol] = castState(eqSS0,eqSS1,p0,p1);
eqTransINBase       = solveTransition('India',p1,eqSS1,initState,initSol,alg);
clear eqSS0 eqSS1 pS0 pS1
%----------------------------------------------------------------------------

disp('Solving US Transition')
initAlg('US');
pUS.qUS2qIndInit = 1.0; % not important
[eqSS0,p0]   = solveBGP(pUS,'US');
[eqSS1,p1]   = solveBGP(pUS,'US');

[initState,initSol] = castState(eqSS0,eqSS1,p0,p1);
eqTransUSBase       = solveTransition('US',p1,eqSS1,initState,initSol,alg);
clear eqSS0 eqSS1 p0 p1

disp('----------------------------------------------')

% -------------------------------------------------------------------------
% --------------------------- COUNTERFACTUAL ------------------------------
% -------------------------------------------------------------------------

disp('COUNTERFACTUAL EXERCISES...')
disp('India: Change  alpha to the US value')
initAlg('India');
pINCF       = pIN; 
pINCF.alpha = pUS.alpha;

% at s.s.
eqINCF           = solveBGP(pINCF,'India');
[~,mINCF,eqINCF] = callMomentBGP(pINCF,eqINCF,1);

disp('Solving India Transition for counterfactual (this will take a while...)')
[initState,initSol] = castState(eqIN,eqINCF,pIN,pINCF);
eqTransINCF   	    = solveTransition('India',pINCF,eqINCF,initState,initSol,alg);
%---------------------------------------------------------------------------

disp('US: Change alpha to Indian level')
initAlg('US');
pUSCF       = pUS; 
pUSCF.alpha = pIN.alpha;

% at s.s.
eqUSCF           = solveBGP(pUSCF,'US');
[~,mUSCF,eqUSCF] = callMomentBGP(pUSCF,eqUSCF,1);
%---------------------------------------------------------------------------

% Other Counterfactuals

disp('Giving Indian delta to US')
p       = pUS;
p.delta = pIN.delta;
eqUSwINdelta                 = solveBGP(p,'US');
[~,mUSwINdelta,eqUSwINdelta] = callMomentBGP(p,eqUSwINdelta,1);
%--------------------------------------------------------------

disp('Giving US delta to India')
p       = pIN;
p.delta = pUS.delta;
eqINwUSdelta                 = solveBGP(p,'India');
[~,mINwUSdelta,eqINwUSdelta] = callMomentBGP(p,eqINwUSdelta,1);
%--------------------------------------------------------------

disp('US economy with other values of delta')
p       = pUS;
p.delta = .75;
eqUSwdelta75                  =  solveBGP(p,'US');
[~,mUSwdelta75,eqUSwdelta75]  =  callMomentBGP(p,eqUSwdelta75,1);

p       = pUS;
p.delta = .999;
eqUSwdelta1                =  solveBGP(p,'US');
[~,mUSwdelta1,eqUSwdelta1] =  callMomentBGP(p,eqUSwdelta1,1);
%------------------------------------------------------------------------

disp('US economy with other values of theta')
p       = pUS;
p.theta = pUS.theta/2;
eqUSwUStheta_05                    =  solveBGP(p,'US');
[~,mUSwUStheta_05,eqUSwUStheta_05] =  callMomentBGP(p,eqUSwUStheta_05,1);

alg.delt = 0.005;
p       = pUS;
p.theta = pUS.theta*2;
eqUSwUStheta_2                   =  solveBGP(p,'US');
[~,mUSwUStheta_2,eqUSwUStheta_2] =  callMomentBGP(p,eqUSwUStheta_2,1);
%------------------------------------------------------------------------

% ADDITIONAL RESULTS FROM BASELINE CALIBRATION
%-----------------------------------------------------------------------
% manager share by alpha for a given firm size distibution (for Figure 5)
initAlg('India')
p = pIN;
alphaGrid = linspace(pIN.alpha, pUS.alpha,2);
resultAlpha = [];
for i = 1:length(alphaGrid)
	p.alpha = alphaGrid(i);
	[res,fval,flag] = fsolve(@(x)solveBGPGivenFSDFunc(x,p,eqIN),[eqIN.omegaM eqIN.omegaP],alg.mopts);

	[~,eqNew] = solveBGPGivenFSDFunc(res,p,eqIN);
	[~,mNew,eqNew] = callMomentBGP(p,eqNew,1);

	resultAlpha = [resultAlpha; p.alpha 1 - eqNew.lProdS [eqNew.manager./eqNew.emp]'];
end
%-----------------------------------------------------------------------

% Product line dist
qq = eqUS.fAll./sum(eqUS.fAll);
productDistUS = [qq(1:9); sum(qq(10:end))];
qq = eqIN.fAll./sum(eqIN.fAll);
productDistIN = [qq(1:9); sum(qq(10:end))];

% Manager share by Size
[managerShareModelBySize,quantiles] = calculateManagerShareBySize(eqIN);

% moment sensitivity
disp('Moment sensitivity')
momSensitivityMatrix = momentSensitivity();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
%--------------------------- ROBUSTNESS ----------------------------------
%-------------------------------------------------------------------------
disp('----------------------------------------------')
disp('----------------------------------------------')
robustness; % runs all the robustness

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------
%--------------------------- TABLES ----------------------------------
%---------------------------------------------------------------------
disp('----------------------------------------------')
disp('TABLES')
table3;
table4;
table5;
table6;
table7;
table8;
table9;
table14;
%---------------------------------------------------------------------
%--------------------------- FIGURES ----------------------------------
%---------------------------------------------------------------------
disp('----------------------------------------------')
disp('FIGURES')
figures;
disp('----------------------------------------------')
disp('DONE! Results are under the folder "output".');

