function [eq]   = ageSizeDistBGP(p,eq,ageMax)

global alg

xStar           = eq.xStar;
tauH            = eq.tauH*ones(alg.NN,1); 
tauL            = eq.tauL;
Tmax            = ageMax/alg.delt;

% Set up some vectors
nvec            = [1:1:alg.NN]';

% Set initial distribution
muHigh          = zeros(alg.NN + 1,1);   % The first row holds exiters
muHigh(2)       = eq.z*p.delta;
distHigh        = zeros(alg.NN+1,Tmax);
distHigh(:,1)   = muHigh;

muLow           = [0,eq.z*(1 - p.delta)]';
distLow         = zeros(2,Tmax);
distLow(:,1)    = muLow;

xStarDelt        = xStar*alg.delt;
tauHDelt         = tauH*alg.delt;
tauLDelt         = tauL*alg.delt;

% Set new probabilities to write in vector notation
u_below         = [0;xStarDelt(1:(alg.NN-1))];
d_above         = [tauHDelt(1:(alg.NN));0];


for s = 1:(Tmax-1)
% future distribution
    mu_primeHigh = muHigh.*(1 - [0;nvec].*([0;xStarDelt] + [0;tauHDelt])) ...           % outflows
                + d_above.*[nvec(1:end);(alg.NN+1)].*[muHigh(2:end);0] ...                      % firms who lose a product
                + [0;u_below].*[0;0;nvec(1:(end-1))].*[0;muHigh(1:(end-1))];                    % firms who win a product
            
    mu_primeLow = muLow.*(1 - [0;(tauLDelt)]) ...
                + [tauLDelt;0].*nvec(1:2).*[muLow(2:end);0]; 
                                          
    distHigh(:,s + 1) = max(mu_primeHigh,0); 
    distLow(:,s + 1)  = max(mu_primeLow,0);
    % update mu for iteration        
    muHigh = max(mu_primeHigh,0);  
    muLow  = max(mu_primeLow,0);
end

% Age Distribution (ageXnumber of prod)
eq.agePoint     = 0:(Tmax-1);
eq.agePoint     = (eq.agePoint*alg.delt)';
eq.distHigh     = distHigh';
eq.distLow      = distLow';
eq.distAll      = eq.distHigh + [eq.distLow zeros(size(eq.distHigh,1),size(eq.distHigh,2) - 2)];


end