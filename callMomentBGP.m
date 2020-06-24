function [score,m,eq] = callMomentBGP(params,eq,calExtra)
  %% Calculates the moment based on equilibrium object

  global alg;
  if ~isstruct(params)
      [p] = putParams(params);    % this part used for estimation
  else
      p  = params;
  end

  score = inf;
  m     = struct();

  % load data moments
  fid = fopen(alg.momentDataFile);
  C   = textscan(fid, '%s%f','delimiter',',');
  fclose(fid);
  m.data = cell2struct(num2cell(C{2}),C{1});

  m.momentModel = zeros(1,length(fieldnames(m.data)));
  m.momentData  = zeros(1,length(fieldnames(m.data)));
  m.momentWgt   = zeros(1,length(fieldnames(m.data)));
  m.momentName  = {};
  mPos          = 1;

  maxAge = alg.maxAgeEst;
  state  = (1:1:alg.NN)';

  % Efficiency units by product
  nStarD    = min(alg.NN(end),ceil(eq.nStar));
  eAll      = zeros(alg.NN,1);

  for j = 1:nStarD-1
      eAll(j) = (p.T/j);
  end

  eAll(nStarD:end) = eq.eM;
  mEff             = (1./(1 - eAll.^p.sigma));

  % Getting grid of number of workers and managers and total employment for n - product firms 
  normWork         = sum(eq.fAll.*(state./mEff))/eq.lProdD;
  worker           = (state./mEff)./normWork;                    % number of worker for firms with n product, it integrates to lprodD over firms.
  managerPerProd   = max((1/p.alpha).*(p.T.*(1/eq.nStar - 1./state)),0);   % manager demand per product for all firms (includes zero demand)
  managerPerFirm   = managerPerProd.*state;
  normManager      = sum(eq.fAll.*managerPerFirm)/(eq.Ltotal - eq.lProdD); 
  manager          = managerPerFirm./normManager;

  if eq.xD ==0
      manager = zeros(alg.NN,1);
  end

  employment        = worker + manager;        % total employment for an n-product firm
  eq.emp            = employment;
  eq.worker         = worker;
  eq.manager        = manager;
  eq.managerPerProd = managerPerProd;
  eq.mEff           = mEff;

  %% Calculate moments

  % 1 - Entry Rate
  m.entRate = eq.z/(eq.lM + eq.hM);  

  % 2 - Exit rates for 1 product firms by age
  exitRateYoungest    = exitRatebyAgeSize1(1,5,1/alg.delt);
  m.exitRatio         = exitRateYoungest/exitRatebyAgeSize1(21,25,1/alg.delt);

  % 3 - Life Cycle 
  meanEmp0_5     = lifeCycle(0,5);
  m.meanEmp21_25 = lifeCycle(21,25)/meanEmp0_5;
  

  % 4 - Share of Employment
  m.shareTot      = aggEmpShare(0,maxAge-1);
  m.shareEmp21_25 = aggEmpShare(21,25)/m.shareTot;

  % 5 - Share of Managers in total labor force
  m.shareMan = max(0,eq.Ltotal - eq.lProdD) ;                

  % 6 - Variance of Log Manager Wage
  m.varLogomegaM     = 1/p.vartheta^2;

  % 7 - Entrepreneur Profit
  m.entProfitShare = eq.fAll'*eq.returnFun;

  % 8 - Share of managerial payment in profit
  max1              = max(0,(1/eq.nStar - 1./state));
  max2              = max(1/eq.nStar, 1./state).^p.sigma;
  m.paymentShareMan = p.sigma*(eq.nStar^(1 - p.sigma))*sum(eq.muAll(1:end).*max1(1:end))/sum(eq.muAll(1:end).*max2(1:end));  
 
  % Moments for calibration
  %--------------------------------------------------------------------------------------------------------------    
  addMoment(m.entRate,        'entRate'        ,alg.wgtentRate);                % moment-1
  addMoment(m.exitRatio,      'exitRatio'      ,alg.wgtexitRatio)               % moment-2    
  addMoment(m.meanEmp21_25,   'meanEmp21_25'   ,alg.wgtmeanEmp21_25);           % moment-3
  addMoment(m.shareEmp21_25,  'shareEmp21_25'  ,alg.wgtshareEmp21_25);          % moment-4
  addMoment(m.shareMan,       'shareMan'       ,alg.wgtshareMan)                % moment-5
  addMoment(m.entProfitShare, 'entProfitShare' ,alg.wgtentProfitShare)          % moment-6
  addMoment(m.varLogomegaM,   'varLogomegaM'   ,alg.wgtvarLogomegaM)            % moment-7
  addMoment(m.paymentShareMan,'paymentShareMan',alg.wgtpaymentShareMan)         % moment-8
  %--------------------------------------------------------------------------------------------------------------    
      

  %% Here we calculate all other moments for reporting
  if calExtra == 1
          ageLb1 = [0  6 11 16 21 26         ];
          ageLb2 = [0  6 11 16 21 26 31 36 41];
          ageLb3 = [0  2  4  6 11 16 21 26 31];

          ageUb1 = [5 10 15 20 25          maxAge-1];
          ageUb2 = [5 10 15 20 25 30 35 40 maxAge-1];
          ageUb3 = [1  3  5 10 15 20 25 30 35];

      
      m.shareEmpOneProdAgg = eq.fAll(1)*employment(1)./eq.Ltotal;
      for j = 1:length(ageLb1)  
          m.meanEmp(j)               = lifeCycle    (ageLb1(j),ageUb1(j));    
          m.meanEmpHigh(j)           = lifeCycleHigh(ageLb1(j),ageUb1(j));
          m.shareOneProd(j)          = shareOneProd(ageLb1(j),ageUb1(j));
          m.shareEmpOneProd(j)       = shareEmpOneProd(ageLb1(j),ageUb1(j));
          m.shareHighAgeBin1(j)      = shareHigh(ageLb1(j),ageUb1(j));
          m.numFirmAgeAgeBin1(j)     = numFirm(ageLb1(j),ageUb1(j)); 

          m.meanProd(j)              = meanProdByAge(ageLb1(j),ageUb1(j));    
          m.meanProdHigh(j)          = meanProdHighByAge(ageLb1(j),ageUb1(j)); 
      end
      
      m.meanEmpNorm 	   = m.meanEmp    ./m.meanEmp(1);
      m.meanEmpNormHigh  = m.meanEmpHigh./m.meanEmpHigh(1);
      
      m.shareEmpHighTot  = max(0,eq.fHigh'*employment)./eq.Ltotal;   
      m.shareOneProdNorm = m.shareOneProd./m.shareOneProd(1);

      for j = 1:length(ageLb3)
          m.shareExitin1ByAgeBySizeOneProdAlt(j) = exitRatebyAgeSize1(ageLb3(j),ageUb3(j),1/alg.delt);
      end

      for j = 0:40
          m.numFirmAge(j+1)           = numFirm(j,j); 
          m.empFirmAge(j+1)           = empFirm(j,j);
          m.shareHigh(j+1)            = shareHigh(j,j);
          m.shareEmpHigh(j+1)         = shareEmpHigh(j,j);
          m.shareEmpOneProd(j+1)      = shareEmpOneProd(j,j);
      end 

      % Aggregate employment in firms without an outside manager 
      nStarD    = min(alg.NN(end),ceil(eq.nStar));
      m.shareEmpAllSE     = max(0,eq.Ltotal - eq.fAll(nStarD:end)'*employment(nStarD:end));   % employment (including managers) share of no outside manager firms 

  end


  %%%%%%%%%%%%%%%
  %% SET SCORE %%
  %%%%%%%%%%%%%%%
  sel            =(m.momentWgt>0);

  m.momentModel  = m.momentModel(sel);
  m.momentData   = m.momentData(sel);
  m.momentWgt    = m.momentWgt(sel);
  m.momentName   = m.momentName(sel);


  m.nMoment = length(m.momentModel);
  m.mErr    = m.momentData - m.momentModel;

  score     = (sum((m.momentWgt'.*(abs(m.mErr')./(0.5*abs(m.momentModel') + 0.5*abs(m.momentData')))).^1))/m.nMoment;

  m.score   = score;
  eq = orderfields(eq);
  m  = orderfields(m); 

  %----------------------------------------------------------------
  %% Functions
  function q = meanProdByAge(minAge,maxAge)
      normDist = eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)./sum(sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)));
      q        = sum(normDist*state);
  end

  function q = meanProdHighByAge(minAge,maxAge)
      normDist = eq.distHigh((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)./sum(sum(eq.distHigh((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)));
      q        = sum(normDist*state);
  end

  function q = lifeCycle(minAge,maxAge)
      normDist = eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)./sum(sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)));
      q        = sum(normDist*employment);
  end

  function q = lifeCycleHigh(minAge,maxAge)
      normDist = eq.distHigh((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)./sum(sum(eq.distHigh((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)));
      q        = sum(normDist*employment);
  end

  function q = numFirm(minAge,maxAge)
      q = sum(sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)));
  end

  function q = empFirm(minAge,maxAge)
      q = sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)*employment);    
  end


  function q = shareHigh(minAge,maxAge)
      q = sum(sum(eq.distHigh((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)))/sum(sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)));
  end

  function q = shareEmpHigh(minAge,maxAge)
      q = sum(sum(eq.distHigh((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)*employment))/sum(sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)*employment));    
  end

  function q = shareOneProd(minAge,maxAge)
      q = sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2))/sum(sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)));
  end

  function q = shareEmpOneProd(minAge,maxAge)
      q = sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2)*employment(1))/sum(sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)*employment));    
  end



  function q = aggEmpShare(minAge,maxAge)
      q = sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2:end)*employment)*alg.delt;
  end


  function q = exitRatebyAgeSize1(minAge,maxAge,N) 
      probH = eq.distHigh(N + 1,1)/eq.distHigh(1,2);
      probL = eq.distLow(N + 1,1)/eq.distLow(1,2);
      
      exitH = sum(eq.distHigh(((1 + (minAge)/alg.delt):(maxAge+1)/alg.delt) + 1,2))*probH;
      exitL = sum(eq.distLow(((1 + (minAge)/alg.delt):(maxAge+1)/alg.delt) + 1,2))*probL;
      
      q = (exitH + exitL)/sum(eq.distAll((1 + minAge/alg.delt):(maxAge+1)/alg.delt,2))/(N*alg.delt);

  end

  %------------------------------------%
  function addMoment(val,str,wgt)
    m.momentModel(mPos) = val;
    m.momentData(mPos)  = m.data.(str); 
    m.momentName{mPos}  = str;
    m.momentWgt(mPos)   = wgt;
    mPos                = mPos + 1;
  end

end
