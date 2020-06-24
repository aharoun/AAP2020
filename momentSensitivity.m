function out = momentSensitivity()

	%  Calculates the percentage change in the moment (row) for a 1% change in the parameter (column) 
	%  from its baseline value, while keeping the rest of the parameters at their benchmark values. 
	%  It report the average of the +1% and -1% changes.

	global alg

	[eqUS,mUS,pUS]  = solver('US');
	[eqIN,mIN,pIN]  = solver('India');
	[meanReg,~,~,~] = bloomExercise(eqIN,pIN,eqUS);

	list1         = [1:2 4:7];
        list2         = [1 3 4 2 5 6];

	[~,pNames]       = readParam(alg.paramsFile);
	mBase 		 = [mIN.momentModel(list2)';meanReg];

	mPlusOnePerc  = [];
	mMinusOnePerc = [];
	pINChange     = pIN;
	pUSChange     = pUS;

	per           = 0.01;


	% PARAMETERS THAT ARE CALIBRATED TO INDIAN 

	for j = 1:length(list1)
	  	
		pINChange.(char(pNames(list1(j)))) = pINChange.(char(pNames(list1(j))))*(1.0 + per);
		[eqINChange,~,~]  =  solveBGP(pINChange,'India');
		[~,mINChange,~]   =  callMomentBGP(pINChange,eqINChange,1);
		
		[meanReg,~,~,~] = bloomExercise(eqINChange,pINChange,eqUS);

		mPlusOnePerc = [mPlusOnePerc [mINChange.momentModel(list2)';meanReg]];
		pINChange   = pIN;

	end

	for j = 1:length(list1)
	  	
		pINChange.(char(pNames(list1(j)))) = pINChange.(char(pNames(list1(j))))*(1.0 - per);
		[eqINChange,~,~]  =  solveBGP(pINChange,'India');
		[~,mINChange,~]   =  callMomentBGP(pINChange,eqINChange,1);
		
		[meanReg,~,~,~] = bloomExercise(eqINChange,pINChange,eqUS);
		
		mMinusOnePerc = [mMinusOnePerc [mINChange.momentModel(list2)';meanReg]];
		pINChange    = pIN;

	end

	perChangePlus     = (mPlusOnePerc  - repmat(mBase,1,length(list1)))./( per*repmat(mBase,1,length(list1)));
	perChangeMinus    = (mMinusOnePerc - repmat(mBase,1,length(list1)))./(-per*repmat(mBase,1,length(list1)));
	perChangeAverage  = (perChangePlus + perChangeMinus)./(2);

	%---------------------------------------------------------------------------------------------------------------
	%---------------------------------------------------------------------------------------------------------------

	% DO sigma separately as it is estimated jointly

	pINChange.sigma = pINChange.sigma*(1.0 + per);
	[eqINChange,~,~]  =  solveBGP(pINChange,'India');
	[~,mINChange,~]   =  callMomentBGP(pINChange,eqINChange,1);

	pUSChange.sigma = pUSChange.sigma*(1.0 + per);
	[eqUSChange,~,~]  =  solveBGP(pUSChange,'US');
	[~,mUSChange,~]   =  callMomentBGP(pUSChange,eqUSChange,1);

	[meanReg,~,~,~] = bloomExercise(eqINChange,pINChange,eqUSChange);

	mPlusOnePercSigma = [mINChange.momentModel(list2)';meanReg];
	pINChange   = pIN;
	pUSChange   = pUS;


	pINChange.sigma = pINChange.sigma*(1.0 - per);
	[eqINChange,~,~]  =  solveBGP(pINChange,'India');
	[~,mINChange,~]   =  callMomentBGP(pINChange,eqINChange,1);

	pUSChange.sigma = pUSChange.sigma*(1.0 - per);
	[eqUSChange,~,~]  =  solveBGP(pUSChange,'US');
	[~,mUSChange,~]   =  callMomentBGP(pUSChange,eqUSChange,1);

	[meanReg,~,~,~] = bloomExercise(eqINChange,pINChange,eqUSChange);

	mMinusOnePercSigma = [mINChange.momentModel(list2)';meanReg];

	perChangePlusSigma     = (mPlusOnePercSigma  - mBase)./( per*mBase);
	perChangeMinusSigma    = (mMinusOnePercSigma - mBase)./(-per*mBase);

	perChangeAverageSigma  = (perChangePlusSigma + perChangeMinusSigma)./(2);


	%---------------------------------------------------------------------------------------------------------------
	%---------------------------------------------------------------------------------------------------------------

	out = [{''} [pNames(list1) 'sigma'];[ [mIN.momentName(list2)';'blooomMoment']  num2cell([perChangeAverage perChangeAverageSigma])]];

end
