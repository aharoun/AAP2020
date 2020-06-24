function [meanReg,meanProfitT,meanProfitNT,flagV,eqINTreat,pINTreat] = bloomExercise(eqIN,pIN,eqUS)
	%% Bloom exercise
	global alg

	eRel   = eqUS.eM/eqIN.eM;
	varphi = log(eRel)/log(1/0.25); 
	lambda = (0.63/0.25)^varphi;
	lambdaTilde = ((lambda*eqIN.eM)^(1-pIN.sigma))*(eqIN.omegaM/(pIN.sigma*pIN.alpha));
	
	
	pINTreat	= pIN;
	pINTreat.alpha  = lambdaTilde*pIN.alpha;
	pINTreat.T      = lambdaTilde*pIN.T;
	
	eqINTreat.tauL   = eqIN.tauL;
	eqINTreat.omegaM = eqIN.omegaM;
	eqINTreat.omegaP = eqIN.omegaP;
	eqINTreat.tauH   = eqIN.tauH;

	[eqINTreat,flagV] = valueFuncBloom(pINTreat,eqINTreat);

	if flagV == 0
	    meanReg = inf;
	    meanProfitT  = 0;
	    meanProfitNT = 0;
	else

	    nYear	  = 107*alg.dtSim;
	    nTreatFirm    = 50;
	    nNonTreatFirm = 50;
	    timeTreat     = 20*alg.dtSim;
	    nCut          = 7;

	    period        = (nYear+timeTreat)/alg.dtSim+1;

	    JMAX   	  = 250;
	    meanProfitT   = zeros(JMAX,period);
	    meanProfitNT  = meanProfitT;
	    regOut        = zeros(JMAX,1);

	    dtSim = alg.dtSim;
	    % firm simulation

	    parfor j = 1:JMAX
		[regOut(j),meanProfitT(j,:),meanProfitNT(j,:)] = bloomRegression(eqIN,eqINTreat,round(j*10000),nYear,...
										nTreatFirm,nNonTreatFirm,timeTreat,nCut,dtSim);
	    end

	    meanReg = mean(regOut);
	end
end
	
	
function [regOut,meanProfitT,meanProfitNT] = bloomRegression(eqIN,eqINTreat,seed,nYear,nTreatFirm,nNonTreatFirm,timeTreat,nCut,dtSim)


	[nProd,exit,treatD] = firmSimBloom(eqIN,eqINTreat,seed,nYear,nTreatFirm,nNonTreatFirm,nCut,dtSim);

	fexitNever  = ~exit;
	indTreat    = find(treatD&fexitNever);
	indNonTreat = find(~treatD&fexitNever);

	indAll      = [indTreat;indNonTreat];
	
	profitTreat    = [repmat(eqIN.piN(nProd(indTreat   ,1)),1,timeTreat/dtSim) eqINTreat.piN(nProd(indTreat,:))];
	profitNonTreat = [repmat(eqIN.piN(nProd(indNonTreat,1)),1,timeTreat/dtSim)      eqIN.piN(nProd(indNonTreat,:))];

	profit = [profitTreat;profitNonTreat];

	profitReg = profit';
	profitReg = profitReg(:);
	
	treatStart = timeTreat/dtSim; 
	treatReg = zeros(length(indAll),nYear/dtSim+treatStart + 1);
	treatReg(1:length(indTreat),treatStart+1:end) = 1;
	treatReg = treatReg';
	treatReg = treatReg(:);

    % Regressions
	logProfit = log(profitReg);
	cons      = ones(length(indAll)*(nYear/dtSim+treatStart +1),1);

	X = 	[treatReg cons];
	
	regOut = (X'*X)\(X'*logProfit);
	% results
	regOut    = regOut(1);
	
	meanProfitT  = mean(profitTreat);
	meanProfitNT = mean(profitNonTreat);	
	
end



