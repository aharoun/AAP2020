function [initState,initSol] = castState(eqSS0,eqSS1,p0,p1)
	%% Initial conditions for transitional dynamics
	global alg
	clear initState
	clear initSol
	
	factor  = inf;

	initState.fH 	    = eqSS0.fHigh;
	initState.fL 	    = 1 - eqSS0.chi;
	initState.qUS2qInd 	= p0.qUS2qIndInit;
	initState.alphaPath = (p1.alpha - p0.alpha)*(1 - exp(-factor*(linspace(0.01,100,alg.T))')) + p0.alpha;

	%% Casting the initial solution to the transition problem

	initSol.r_gY      = (eqSS1.r_gY   - eqSS0.r_gY)  *(1 - exp(-factor*(linspace(0.01,100,alg.T))')) + eqSS0.r_gY;
	initSol.tauH	  = (eqSS1.tauH   - eqSS0.tauH)  *(1 - exp(-factor*(linspace(0.01,100,alg.T))')) + eqSS0.tauH;    
	initSol.omegaM 	  = (eqSS1.omegaM - eqSS0.omegaM)*(1 - exp(-factor*(linspace(0.01,100,alg.T))')) + eqSS0.omegaM;
	initSol.omegaP    = (eqSS1.omegaP - eqSS0.omegaP)*(1 - exp(-factor*(linspace(0.01,100,alg.T))')) + eqSS0.omegaP;

end
