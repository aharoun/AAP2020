function [eq,flagV] = valueFuncBGP(p,eq)
%% Calculates Vn 
% September 2016
global alg


nEnd    = alg.NN;          % maximum number of products
maxErr  = 1.0e-8;
maxIter = 1000;


state  = (1:nEnd)';
  
nStarD          = min(floor(eq.nStar),alg.NN);
piN             = p.T*p.sigma*((eq.nStar/p.T).^(1 - p.sigma)) + (1 - p.sigma)*((p.T/eq.nStar)^p.sigma).*state;
piN(1:nStarD)   = (p.T^p.sigma)*(([1:nStarD]').^(1 - p.sigma));
VnOld           = piN/(p.rho + eq.tauH);     

VnNew           = zeros(alg.NN,1);
V0              = 0;   
    
delta = 100;    

% Now the value function iteration
crit = 1;
iter = 0;
flagV= 1; 

while crit > maxErr && iter<=maxIter  

    xn        = max(zeros(nEnd,1),(p.zeta*(1/eq.omegaP)*(p.theta^(1/p.zeta)).*([VnOld(2:end);2*VnOld(end)-VnOld(end-1) ] - VnOld)).^(p.zeta/(1 - p.zeta))); % xn is calculated given VnOld
    returnFun = max(0,(piN - eq.omegaP*(p.theta^(-1/p.zeta))*(xn.^(1/p.zeta)).*state));    % return function is calculated given VnOld 

    % Implicit method    
    AAp = spdiags( [[state(2:end);1].*(eq.tauH)  -state.*(eq.tauH+xn) [1;state(1:end-1)].*[1;xn(1:end-1)]],-1:1,alg.NN,alg.NN);
    
    AAp(end,end)   = nEnd*(xn(end) - eq.tauH);
    AAp(end,end-1) = -nEnd*(xn(end) - eq.tauH);

    BBp = (1/delta + p.rho)*speye(alg.NN) - AAp;
    bbp = (1/delta)*VnOld + returnFun;


    VnNew = BBp\bbp;
    
    crit  = max((VnNew - VnOld).^2);
    if isnan(crit)
        crit = 1;
    end
    iter = iter + 1;
 
    VnOld = VnNew;
 
end

    if iter < maxIter      
        eq.xStar  = xn;
        eq.VnH    = VnNew;
        eq.VnL    = piN(1)/(p.rho + eq.tauL);
        eq.returnFun = returnFun;
        eq.piN = piN;
        flagV     = 1;
    else
        flagV     = 0;
    end

end
