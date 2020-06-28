function [nProd,exit,treatD] = firmSimBloom(eq,eqTreat,seed,nYear,nTreatFirm,nNonTreatFirm,nCut,dtSim)

rng(seed,'simdTwister');       


Tmax  = nYear/dtSim;

nFirm = nTreatFirm + nNonTreatFirm; 

% Preallocation for reporting
nProd    = zeros(nFirm,Tmax + 1);
exit     = false(nFirm,1);

% Initialize at stationary dist.
nProdD   = zeros(nFirm,1);
treatD   = false(nFirm,1);
treatD(1:nTreatFirm) = true;


nProdD     = nCut;               

nProd(:,1)  = nProdD; 

xDelt       = eq.xStar*dtSim;
xTreatDelt  = eqTreat.xStar*dtSim;
tauHDelt    = eq.tauH*dtSim;


for t = 1:Tmax
    draw          = rand(nFirm,1); 
    stateUp       = ((1 - treatD).*xDelt(nProdD) + treatD.*xTreatDelt(nProdD)).*nProdD;      
    stateDown     = stateUp + tauHDelt.*nProdD ;                  % creative destruction rate is same for all firms, products.
    
    up            = (draw < stateUp);            
    down          = logical((1 - up).*(draw < stateDown));
    nProdD        = nProdD + up - down;                         % update number of product line

    % Check if exit
    ifExit           = (nProdD<1);
    nProdD(ifExit)   = 1;        
    exit(ifExit)     = true;
    
    % reporting
    nProd(:,t+1)= nProdD; 

end
        





end


