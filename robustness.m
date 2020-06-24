%% Robustness exercises
%----------------------------------------------------------------------------------------------------
disp('ROBUSTNESS...')
%----------------------------------------------------------------------------------------------------
%% 1- Estimation with Managerial Payment
   [eqUSWPayment,mUSWPayment,pUSWPayment] = solver('USWPayment');
   [eqINWPayment,mINWPayment,pINWPayment] = solver('IndiaWPayment');

   lambdaMethodWPayment = calibrateDiffusion(pINWPayment,eqINWPayment,pUSWPayment,eqUSWPayment);
   pINWPayment.lambda   = lambdaMethodWPayment(1);

   pINWPaymentCF       = pINWPayment;
   pINWPaymentCF.alpha = pUSWPayment.alpha;

   eqINWPaymentCF                   = solveBGP(pINWPaymentCF,'IndiaWPayment');
   [~,mINWPaymentCF,eqINWPaymentCF] = callMomentBGP(pINWPaymentCF,eqINWPaymentCF,1);


%% 2- Estimation with Bloom and Managerial Payment
   [eqUSWBloomAndPayment,mUSWBloomAndPayment,pUSWBloomAndPayment] = solver('USWBloomAndPayment');
   [eqINWBloomAndPayment,mINWBloomAndPayment,pINWBloomAndPayment] = solver('IndiaWBloomAndPayment');

   pINWBloomAndPaymentCF       = pINWBloomAndPayment;
   pINWBloomAndPaymentCF.alpha = pUSWBloomAndPayment.alpha;

   eqINWBloomAndPaymentCF                           = solveBGP(pINWBloomAndPaymentCF,'IndiaWBloomAndPayment');
   [~,mINWBloomAndPaymentCF,eqINWBloomAndPaymentCF] = callMomentBGP(pINWBloomAndPaymentCF,eqINWBloomAndPaymentCF,1);

%% 3 - Elastic Labor Supply
% with elastic labor supply on a grid
    alg.elasticLaborSupply = 1;

    pINels102  = pIN;
    eqINels102 = eqIN;
    mINels102  = mIN;

    pINels102CF = pINels102;
    pINels102CF.alpha = pUS.alpha;
    pINels102CF.Ltotal = 1.02;

    eqINels102CF                 = solveBGP(pINels102CF,'India');
    [~,mINels102CF,eqINels102CF] = callMomentBGP(pINels102CF,eqINels102CF,1);

    eqUSels102 = eqUS;
    mUSels102 = mUS;

    %---------------------------------------------------------------------------

    pINels105 = pIN;
    eqINels105 = eqIN;
    mINels105 = mIN;

    pINels105CF = pINels105;
    pINels105CF.alpha = pUS.alpha;
    pINels105CF.Ltotal = 1.05;

    eqINels105CF         = solveBGP(pINels105CF,'India');
    [~,mINels105CF,eqINels105CF] = callMomentBGP(pINels105CF,eqINels105CF,1);

    eqUSels105 = eqUS;
    mUSels105 = mUS;

    %---------------------------------------------------------------------------

    alg.elasticLaborSupply = 0;

%% 3 - Entry elasticity
% Low
    pINzetaEntLow            = pIN;
    pINzetaEntLow.zetaEnt    = 0.4;
    pINzetaEntLow.thetaEnt   = getThetaEnt(pINzetaEntLow,eqIN);    % recalibrate thetaEnt s.t. it yields same equilibrium with baseline

    eqINzetaEntLow                   = solveBGP(pINzetaEntLow,'India');
    [~,mINzetaEntLow,eqINzetaEntLow] = callMomentBGP(pINzetaEntLow,eqINzetaEntLow,1);

    pINzetaEntLowCF       = pINzetaEntLow;
    pINzetaEntLowCF.alpha = pUS.alpha;

    eqINzetaEntLowCF                     = solveBGP(pINzetaEntLowCF,'India');
    [~,mINzetaEntLowCF,eqINzetaEntLowCF] = callMomentBGP(pINzetaEntLowCF,eqINzetaEntLowCF,1);

    eqUSzetaEntLow = eqUS;  % just reporting to work below
    mUSzetaEntLow  = mUS;
% High
    pINzetaEntHigh             = pIN;
    pINzetaEntHigh.zetaEnt     = 0.6;
    pINzetaEntHigh.thetaEnt    = getThetaEnt(pINzetaEntHigh,eqIN);    % recalibrate thetaEnt s.t. it yields same equilibrium with baseline

    eqINzetaEntHigh                    = solveBGP(pINzetaEntHigh,'India');
    [~,mINzetaEntHigh,eqINzetaEntHigh] = callMomentBGP(pINzetaEntHigh,eqINzetaEntHigh,1);

    pINzetaEntHighCF         = pINzetaEntHigh;
    pINzetaEntHighCF.alpha = pUS.alpha;

    eqINzetaEntHighCF                      = solveBGP(pINzetaEntHighCF,'India');
    [~,mINzetaEntHighCF,eqINzetaEntHighCF] = callMomentBGP(pINzetaEntHighCF,eqINzetaEntHighCF,1);
    
    eqUSzetaEntHigh = eqUS;  % just reporting to work below
    mUSzetaEntHigh  = mUS;
%----------------------------------------------------------------------------------------------------
%% 4 - Innovation elasticity (re-estimated)
% zeta = 0.4
    [eqUSzeta04,mUSzeta04,pUSzeta04] = solver('USzeta04');
    [eqINzeta04,mINzeta04,pINzeta04] = solver('Indiazeta04');

    pINzeta04CF       = pINzeta04;
    pINzeta04CF.alpha = pUSzeta04.alpha;

    eqINzeta04CF                 = solveBGP(pINzeta04CF,'Indiazeta04');
    [~,mINzeta04CF,eqINzeta04CF] = callMomentBGP(pINzeta04CF,eqINzeta04CF,1);
% zeta = 0.6
    [eqUSzeta06,mUSzeta06,pUSzeta06] = solver('USzeta06');
    [eqINzeta06,mINzeta06,pINzeta06] = solver('Indiazeta06');

    pINzeta06CF             = pINzeta06;
    pINzeta06CF.alpha     = pUSzeta06.alpha;

    eqINzeta06CF                 = solveBGP(pINzeta06CF,'Indiazeta06');
    [~,mINzeta06CF,eqINzeta06CF] = callMomentBGP(pINzeta06CF,eqINzeta06CF,1);
%----------------------------------------------------------------------------------------------------

    [eqUSFirm,mUSFirm,pUSFirm] = solver('USFirm');
    [eqINFirm,mINFirm,pINFirm] = solver('IndiaFirm');

    pINFirmCF             = pINFirm;
    pINFirmCF.alpha     = pUSFirm.alpha;

    eqINFirmCF               = solveBGP(pINFirmCF,'IndiaFirm');
    [~,mINFirmCF,eqINFirmCF] = callMomentBGP(pINFirmCF,eqINFirmCF,1);
%----------------------------------------------------------------------------------------------------
%% 7 - Transition speed
% We vary speed of transition parameter 'lambda' such that half life of qIN2qUS is 25% shorter and longer

    halfLifeBase = halfLife(eqIN.tau,eqUS.g,pIN.lambda,pIN.gamma,qUS2qIndt0);
% half-life 25% shorter
    func         = @(x)halfLife(eqIN.tau,eqUS.g,x,pIN.gamma,qUS2qIndt0) - halfLifeBase*.75;  
    lambdaShort     = fsolve(func,pIN.lambda,alg.mopts); 

    pINHLShort     = pIN;
    pINHLShort.lambda = lambdaShort;

    eqINHLShort                = solveBGP(pINHLShort,'India');
    [~,mINHLShort,eqINHLShort] = callMomentBGP(pINHLShort,eqINHLShort,1);

    pINHLShortCF         = pINHLShort;
    pINHLShortCF.alpha = pUS.alpha;

    eqINHLShortCF                  = solveBGP(pINHLShortCF,'India');
    [~,mINHLShortCF,eqINHLShortCF] = callMomentBGP(pINHLShortCF,eqINHLShortCF,1);

    eqUSHLShort = eqUS;  % just reporting to work below
% half-life 25% longer
    func         = @(x)halfLife(eqIN.tau,eqUS.g,x,pIN.gamma,qUS2qIndt0) - halfLifeBase*1.25;  
    lambdaLong   = fsolve(func,pIN.lambda,alg.mopts);

    pINHLLong     = pIN;
    pINHLLong.lambda = lambdaLong;

    eqINHLLong                = solveBGP(pINHLLong,'India');
    [~,mINHLLong,eqINHLLong]  = callMomentBGP(pINHLLong,eqINHLLong,1);

    pINHLLongCF       = pINHLLong;
    pINHLLongCF.alpha = pUS.alpha;

    eqINHLLongCF                 = solveBGP(pINHLLongCF,'India');
    [~,mINHLLongCF,eqINHLLongCF] = callMomentBGP(pINHLLongCF,eqINHLLongCF,1);

    eqUSHLLong = eqUS;  % just reporting to work below
%----------------------------------------------------------------------------------------------------
%% 8 - High vartheta

    [eqUSvarthetaHigh,mUSvarthetaHigh,pUSvarthetaHigh] = solver('USvarthetaHigh');
    [eqINvarthetaHigh,mINvarthetaHigh,pINvarthetaHigh] = solver('IndiavarthetaHigh');

    pINvarthetaHighCF       = pINvarthetaHigh;
    pINvarthetaHighCF.alpha = pUSvarthetaHigh.alpha;

    eqINvarthetaHighCF                      = solveBGP(pINvarthetaHighCF,'IndiavarthetaHigh');
    [~,mINvarthetaHighCF,eqINvarthetaHighCF] = callMomentBGP(pINvarthetaHighCF,eqINvarthetaHighCF,1);


%% Reporting
    robustnessList =  {'','WPayment','WBloomAndPayment','zetaEntLow','zetaEntHigh','zeta04','zeta06','Firm','HLLong','HLShort','els102','els105','varthetaHigh'};

    rb.tauIN   = zeros(numel(robustnessList),1);
    rb.tauINCF = rb.tauIN;
    rb.mEffIN  = rb.tauIN;
    rb.mEffINCF= rb.tauIN;
    rb.lPIN    = rb.tauIN;
    rb.lPINCF  = rb.tauIN;
    rb.tauUS   = rb.tauIN;
    rb.mEffUS  = rb.tauIN;
    rb.lPUS    = rb.tauIN;
    rb.qIN2US  = rb.tauIN;
    rb.qIN2USCF= rb.tauIN;
    rb.mEmp2125= rb.tauIN;
    rb.mEmp2125CF        = rb.tauIN;
    rb.shareOneProd2125  = rb.tauIN;
    rb.shareOneProd2125CF= rb.tauIN;


    for i = 1:numel(robustnessList)
        rb.tauIN(i)     = eval(['eqIN' robustnessList{i}   '.tau']);
        rb.tauINCF(i)   = eval(['eqIN' robustnessList{i} 'CF.tau']);
        rb.mEffIN(i)    = eval(['eqIN' robustnessList{i}   '.M']);      
        rb.mEffINCF(i)  = eval(['eqIN' robustnessList{i} 'CF.M']);      
        rb.lPIN(i)      = eval(['eqIN' robustnessList{i}   '.lProdD']);   
        rb.lPINCF(i)    = eval(['eqIN' robustnessList{i} 'CF.lProdD']);

        rb.tauUS(i)     = eval(['eqUS' robustnessList{i}   '.tau']);
        rb.mEffUS(i)    = eval(['eqUS' robustnessList{i}   '.M']);
        rb.lPUS(i)      = eval(['eqUS' robustnessList{i}   '.lProdD']);

        parEps       = eval(['pIN' robustnessList{i} '.lambda'   ]);
        parLambda    = eval(['pIN' robustnessList{i} '.gamma']);
        rb.qIN2US(i)    = 1/exp((log(parLambda)/parEps)*(rb.tauUS(i) - rb.tauIN(i))  /rb.tauIN(i));
        rb.qIN2USCF(i)  = 1/exp((log(parLambda)/parEps)*(rb.tauUS(i) - rb.tauINCF(i))/rb.tauINCF(i));

        lM   = eval(['eqIN' robustnessList{i}   '.lM']); 
        hM   = eval(['eqIN' robustnessList{i}   '.hM']); 
        lMCF = eval(['eqIN' robustnessList{i}   'CF.lM']); 
        hMCF = eval(['eqIN' robustnessList{i}   'CF.hM']); 
        rb.aveFirmSize(i)  = 1/(lM +hM);
        rb.aveFirmSizeCF(i)= 1/(lMCF + hMCF);

        rb.shareOneProd2125(i)  = eval(['mIN' robustnessList{i}   '.shareOneProd(5)']);
        rb.shareOneProd2125CF(i)= eval(['mIN' robustnessList{i} 'CF.shareOneProd(5)']);
    end

    rb.yIN2US   = (rb.mEffIN  ./rb.mEffUS).*(rb.lPIN  ./rb.lPUS).*rb.qIN2US;
    rb.yIN2USCF = (rb.mEffINCF./rb.mEffUS).*(rb.lPINCF./rb.lPUS).*rb.qIN2USCF;

    rb.dTauIN            = (rb.tauINCF./rb.tauIN-1)*100;
    rb.dQIN2US           = (rb.qIN2USCF./rb.qIN2US-1)*100;
    rb.dYIN2US           = (rb.yIN2USCF./rb.yIN2US-1)*100;
    rb.dAverageFirmSize  = (rb.aveFirmSizeCF'./rb.aveFirmSize'-1)*100;
    rb.dShareOneProd2125 = (rb.shareOneProd2125CF./rb.shareOneProd2125-1)*100;
