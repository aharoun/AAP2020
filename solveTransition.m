function eqTrans  = solveTransition(cnt,params,eqSS,initState,initSol,alg)
    % solves the transition
    if ~isstruct(params)
        p = putParams(params);
    else
        p = params;
    end

    maxIter      = 10000;
    convergCrit  = 1.0e-5;
    relax        = .9225;   				% for convergence of r-gY and tau
    state        = [1:alg.NN]';
    R            = zeros(alg.NN,1);
    R(1)         = 1;					% selection matrix for firm size distribution

    %% Initialize solution
    r_gY0   = initSol.r_gY;          
    tauH0   = initSol.tauH;        
    omegaM0 = initSol.omegaM;
    omegaP0 = initSol.omegaP;     

    tauHUpdated = zeros(alg.T,1);

    r_gYNew   = r_gY0;
    tauHNew   = tauH0;
    omegaMNew = omegaM0;
    omegaPNew = omegaP0;



    % preallocation
    vHt       = zeros(alg.NN,alg.T);
    vLt       = zeros(1,alg.T);
    xnt       = zeros(alg.NN,alg.T);
    taut      = zeros(alg.T,1);
    AApt      = cell(alg.T,1);  					% save transition matrix for later
    fHt       = cell(alg.T + 1,1);		
    fHt{1}    = initState.fH;   			        % Initial distribution	
    fLt       = zeros(1,alg.T + 1);      
    fLt(1)    = initState.fL;                     	% Initial distribution  
    fAllt     = cell(alg.T + 1,1);
    fAllt{1}  = fHt{1} + [fLt(1);zeros(alg.NN-1,1)];
    lProdD    = zeros(alg.T,1);
    lProdS    = zeros(alg.T,1);
    RDCost    = zeros(alg.T,1);
    entryCost = zeros(alg.T,1);
    Mt        = zeros(alg.T,1);
    managerS  = zeros(alg.T,1);
    managerD  = managerS; 

    p.alphat = initState.alphaPath;     

    iter = 1;
    crit = 1;

    while crit>convergCrit && iter<maxIter
        r_gYt   = r_gYNew;
        tauHt   = tauHNew;
        omegaM  = omegaMNew;
        omegaP  = omegaPNew; 
        
        tauLt = p.beta*tauHt;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% STEP 1: Solve the value function backward in time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nStar    = min(alg.NN-1,max(0.000001,p.T*(omegaM./(p.sigma.*p.alphat)).^(1/(1 - p.sigma))));
        nStarD   = floor(nStar);
        eM       = p.T./nStar;

    	vHCurrent = eqSS.VnH;
        vLCurrent = eqSS.VnL;

        for t = alg.T:-1:1
        	vHt(:,t) = vHCurrent;
            vLt(t)   = vLCurrent;
            
        	pi_n = p.T*p.sigma*((nStar(t)/p.T).^(1 - p.sigma)) + state.*((1 - p.sigma)*((p.T/nStar(t))^p.sigma));
    		pi_n(1:nStarD(t)) = (p.T^p.sigma)*(([1:nStarD(t)]').^(1 - p.sigma));

    		xn  = max(zeros(alg.NN,1),(p.zeta*(1/omegaP(t))*(p.theta^(1/p.zeta))...
                 .*([vHCurrent(2:end);2*vHCurrent(end)-vHCurrent(end-1) ] - vHCurrent)).^(p.zeta/(1 - p.zeta))); % xnt is calculated given vHCurrent
        	
            returnFun = max(0,(pi_n - omegaP(t)*(p.theta^(-1/p.zeta))*(xn.^(1/p.zeta)).*state));    % net return function 
            
        	AAp 		   = spdiags([[state(2:end);1].*(tauHt(t))  -state.*(tauHt(t)+xn) [1;state(1:end-1)].*[1;xn(1:end-1)]],-1:1,alg.NN,alg.NN);
        	AAp(end,end)   =  alg.NN*(xn(end) - tauHt(t));
        	AAp(end,end-1) = -alg.NN*(xn(end) - tauHt(t));

        	AApt{t} = AAp;

        	BBp = (1/alg.dt + r_gYt(t))*speye(alg.NN) - AAp;
        	bbp = (1/alg.dt)*vHCurrent + returnFun;

        	vHCurrent = BBp\bbp;  % high type
            vLCurrent = (vLCurrent + alg.dt*pi_n(1))/(1 + (r_gYt(t) + tauLt(t))*alg.dt);
            xnt(:,t)  = xn;
    	end
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% STEP 2: Solve the firm size distribution forward in time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % get the entry rate 
        zt = (p.thetaEnt^(1/(1 - p.zetaEnt))).*((1./omegaP).^(p.zetaEnt/(1 - p.zetaEnt)))...
            .*((p.zetaEnt)^(p.zetaEnt/(1 - p.zetaEnt))).*((p.delta*vHt(1,:)' + (1 - p.delta)*vLt').^(p.zetaEnt/(1 - p.zetaEnt)));        

        for t = 1:alg.T
        	AApTT      = AApt{t}';

        	fHtPlus1   = (speye(alg.NN) - AApTT*alg.dt)\(fHt{t} + alg.dt*R*p.delta*zt(t));    		% implicit method
            fLt(t+1)   = (fLt(t) + alg.dt*(1 - p.delta)*zt(t))/(1 + alg.dt*tauLt(t));
            fHt{t+1}   = fHtPlus1;
            fAllt{t+1} = fHtPlus1 + [fLt(t+1);zeros(alg.NN-1,1)]; 

            taut(t)    = sum(state.*xnt(:,t).*fHt{t}) + zt(t);				  					% aggregate tau
        end

        tauHUpdated    = taut./(1 + fLt(1:end-1)'*(p.beta - 1));
    	
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% STEP 3: Check demand and supply for production workers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for t = 1:alg.T
        	% Manager Demand
    		xDAux       = max((1/p.alphat(t)).*(p.T.*(1./nStar(t) - 1./state)).*fAllt{t}.*state,0);
    		xD          = sum(xDAux); 
            managerD(t) = xD;
    		% Aggregate effort
    		nStarF    = min(alg.NN,ceil(nStar(t)));

    		partA     = 0;
    		for i=1:(nStarF - 1)
        		partA = partA + fAllt{t}(i)*i*(1 - (p.T/i)^p.sigma);
    		end
    		partB      = (1 - eM(t)^p.sigma)*sum(fAllt{t}(nStarF:end).*state(nStarF:end));
    		M          = ((partA + partB)^(-1));

    		% Production worker demand and supply
    		lProdS(t)  = max(0,(1 - ((p.muM*(p.vartheta - 1)/p.vartheta)^p.vartheta)*((omegaM(t)/omegaP(t))^(p.vartheta))));
            Mt(t)      = M; 

    		RDCost(t)    = omegaP(t)*sum((p.theta^(-1/p.zeta))*(xnt(:,t).^(1/p.zeta)).*state.*fHt{t});  % RD cost over Y
    		entryCost(t) = omegaP(t)*(p.thetaEnt^(-1/p.zetaEnt))*(zt(t).^(1/p.zetaEnt));			       % Entry cost over Y	
    	end


        % updated wages based on demand=supply condition
        A1            = (p.vartheta/(p.vartheta - 1))*((p.muM*(p.vartheta - 1)/p.vartheta)^p.vartheta);
        omegaMUpdated = (managerD./A1).^(1/(p.vartheta-1)).*omegaP;
        omegaPUpdated = 1./(lProdS.*Mt);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% STEP 4: get the path for step size, stepSize and corresponding quality growth rate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(cnt, 'US')
            stepSize = p.gamma;    
            gQt      = log(stepSize).*taut;
        elseif strcmp(cnt, 'India')
            stepSize = getStepSizePath(taut,eqSS.g,initState.qUS2qInd,p,alg);
            gQt      = log(stepSize(1:end-1)').*taut;
        end
                                     
        

    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% STEP 5: Calculate consumption growth and get the new r-gY difference
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        qPath  = cumprod([1; (1 + gQt(1:end-1)*alg.dt)]);
        yPath  = Mt.*lProdS.*qPath; % this is as if Yt = MtxLtxQt-1
        gYt    = (yPath(2:end) - yPath(1:end-1))./(yPath(1:end-1)*alg.dt);
        gYt    = [gYt;gYt(end)];        % to make it Tx1 size, under the assp. that it reaches ss.
        
        cPath  = yPath.*(1 - RDCost - entryCost);
        gCt    = (cPath(2:end) - cPath(1:end-1))./(cPath(1:end-1)*alg.dt);
        gCt    = [gCt;gCt(end)];  		% to make it Tx1 size, under the assp. that it reaches ss.

        r_gYUpdated = gCt + p.rho - gYt;
        rt          = r_gYt + gYt;

        %%%%%%%%%%%%%%%%%%%%%%%%%
        %% Convergence criteria
        %%%%%%%%%%%%%%%%%%%%%%%%%
        r_gYNew    = relax.*r_gYt  + (1 - relax).*r_gYUpdated;
        tauHNew    = relax.*tauHt  + (1 - relax).*tauHUpdated;
        omegaPNew  = relax.*omegaP + (1 - relax).*omegaPUpdated;
        omegaMNew  = relax.*omegaM + (1 - relax).*omegaMUpdated;


        

        crit     = max(abs([r_gYt - r_gYUpdated;tauHt - tauHUpdated;omegaP - omegaPUpdated;omegaM - omegaMUpdated]));
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mod(iter,10)==0
            fprintf('Iteration = %3i, Criteria = %8.2E\n',iter,crit);
    	end

        critSeq(iter) = crit;
        iter = iter+1;

    end


    disp('-> CONVERGED!')
  
    eqTrans.xt       = xnt;
    eqTrans.zt       = zt;
    eqTrans.omegaM   = omegaM;
    eqTrans.omegaP   = omegaP;
    eqTrans.tau      = taut;
    eqTrans.tauH     = tauHt;
    eqTrans.tauL     = tauLt;
    eqTrans.stepSize = stepSize';
    eqTrans.gY       = gYt;
    eqTrans.gC       = gCt;
    eqTrans.gQ       = gQt;
    eqTrans.lProdS   = lProdS;
    eqTrans.fHt      = fHt;
    eqTrans.fLt      = fLt;
    eqTrans.fAllt    = fAllt;
    eqTrans.nStar    = nStar;
    eqTrans.alphat   = p.alphat;
    eqTrans.Mt       = Mt;
    eqTrans.RDCost   = RDCost;
    eqTrans.entryCost= entryCost;
    eqTrans.r_gY     = r_gYt; 


end
