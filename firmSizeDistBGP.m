function [eq,flag] = firmSizeDistFast(params,eq)

%% Solves stationary distribution given the value of xStar and structural parameters

% f(n)  : share of firms with n product (NOT normalized to 1)
% mu(n) : share of product lines which belongs to a firm with n-product (sum is normalized to one for aggregate)
% lM + hM gives total mass of firms, given that total number of product line is normalized to one.

global alg
if ~isstruct(params)
    p = putParams(params);    % this part used for estimation
else
    p = params;
end


state     = (1:1:alg.NN)';
EL        = eq.z*(1 - p.delta); % entry to low type
EH        = eq.z*p.delta;       % entry to high type 
   
options = optimset('Display','off','MaxFunEvals',2000,'MaxIter',1000,'TolFun', 1e-8, 'TolX', 1e-8);
[result,f,flag] = fsolve(@objFunc, eq.tauH*1.5,options); 
 if flag<0
     %disp('Second solver for Firm Size dist')
     options = optimset('Display','off','MaxFunEvals',2000,'MaxIter',1000,'TolFun', 1e-7, 'TolX', 1e-7);
     objFunc2 = @(x)sum(objFunc(x).^2);
     [result,f,flag] = fminsearch(objFunc2,eq.tauH,options);
     if f>1
         warning('Could not solve the FSD...')
     end
     [result,f,flag] = fsolve(@objFunc, result,options);
 end
 if flag<0
     warning('Could not solve FSD...')
 else

 
[~,eqF]    = objFunc(result);
eq.tau     = result; 
eq.chi     = 1 - eqF.fLow(1);
eq.tauHNew = eq.tau/(eq.chi + (1 - eq.chi)*p.beta);
eq.tauLNew = eq.tauHNew*p.beta;
eq.hM      = sum(eqF.fHigh);
eq.lM      = eqF.fLow(1);
eq.fHigh   = eqF.fHigh;
eq.fAll    = (eqF.fHigh + eqF.fLow);              
eq.muHigh  = eqF.fHigh.*state;
eq.muLow   = eqF.fLow.*state;
eq.muAll   = eq.muHigh + eq.muLow ;



 end
function [diff,eqF,flagT] = objFunc(x)

    tau = x;

    tauH       = tau + EL*((1 - p.beta)/p.beta);
    tauL       = tauH*p.beta;
    xH         = eq.xStar;

    eqF.fHigh   = max(([1;cumprod(xH(1:end-1))].*EH)./(state.*tauH.^state),0);
    eqF.fHigh(isnan(eqF.fHigh)) = 0;
    eqF.fHigh(isinf(eqF.fHigh)) = 0;
    eqF.fLow    = zeros(alg.NN,1);
    eqF.fLow(1) = EL/tauL;

    tauNew     = sum(state.*xH.*eqF.fHigh) + EL + EH;
    
    diff = tau-tauNew;
    
    if isnan(diff) || isinf(diff)
        diff = 100000;
        eqF  = 0; 
    end
    end
    
end
% end of function