function initAlg(cnt)
  global alg
  warning off;

  alg.cnt = 'US';
  if (nargin == 1)
    alg.cnt = cnt;
  end

  %% Parameters for solver
  alg.NN            = 200; 
  alg.maxIter       = 1000;    
  alg.maxIter2      = 100; 
  alg.maxIter3      = 50; 
  alg.crit1         = 1e-6;
  alg.crit2         = 1e-9;
  alg.final         = 0;
  alg.saveInit      = 0;

  alg.mopts = optimset('Display','off','MaxFunEvals',300,'MaxIter',300,'TolX',1.0e-12,'TolFun',1.0e-12);
  
  %% transition
  alg.tYear = 400;
  alg.dt    = .1;
  alg.dtSim = 1/50;
  alg.T     = alg.tYear/alg.dt;
  alg.nPool = 4;
  
  %% Moments
  alg.maxAgeEst     = 100;                  
  alg.maxAgePostEst = 100;                  
  alg.delt          = 1/50;               
    
  %% File names

  if strcmp(alg.cnt(1:2), 'US')
    if any(strfind(alg.cnt,'Firm'))
      momentcnt = 'USFirm';
    else
      momentcnt = 'US';
    end
  elseif strcmp(alg.cnt(1:2), 'In')
    momentcnt = 'India';
  end
    
  alg.paramsFile      = ['params'  filesep 'params' '-' alg.cnt '.txt'];
  alg.eqvarsFile      = ['eqvars'  filesep 'eqvars' '-' alg.cnt '.txt'];
  alg.momentDataFile  = ['moments' filesep 'moments' '-' momentcnt '.csv'];
  
  %% parameters for transition
  alg.opt23       = odeset('RelTol',1e-4 ,'AbsTol',1e-4,'Stats','off','Vectorized','on');
  alg.timePath    = linspace(alg.dt,alg.tYear,alg.T)';
  alg.timePath2   = linspace(alg.dt,alg.tYear,5)';

  alg.elasticLaborSupply = 0;

  %% ESTIMATION WEIGHTS
  alg.wgtmeanEmp21_25   = 1.0;
  alg.wgtshareMan       = 1.0;
  alg.wgtentRate        = 1.0;
  alg.wgtentProfitShare = 1.0;
  alg.wgtshareEmp21_25  = 1.0;
  alg.wgtexitRatio      = 1.0;

  alg.wgtpaymentShareMan= 0.0;
  alg.wgtvarLogomegaM   = 0.0;
      
  if any(strfind(alg.cnt,'Payment'))
    alg.wgtpaymentShareMan =  1.0;                                                 
  end


end
