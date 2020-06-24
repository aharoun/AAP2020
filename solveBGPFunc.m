function [diff,eq] = solveBGPFunc(x,p)
%% Solves equilibrium given parameters
global alg

eq.omegaM     = x(1);
eq.omegaP     = x(2);
eq.tauH       = x(3);

eq.tauL       = p.beta*eq.tauH; 
eq.nStar      = p.T*(eq.omegaM/(p.sigma*p.alpha))^(1/(1 - p.sigma));

if eq.nStar <= p.T || isnan(eq.nStar) || isinf(eq.nStar) || ~isreal(eq.nStar) || eq.omegaP<0 || isnan(eq.omegaP) || isinf(eq.omegaP) || ~isreal(eq.omegaP)
    
    diff = ones(3,1)*10000;
    eq   = 0.0;
else

	state     = (1:1:alg.NN)';

	[eq,flagV] = valueFuncBGP(p,eq);

	if flagV == 0
	    diff = ones(3,1)*100000;
	else
	    	% Total labor supply
		if alg.elasticLaborSupply == 1
		    eq.Ltotal = p.Ltotal;
		else
		    eq.Ltotal =  1.0;   % this is our benchmark model
		end
	    
		%% Step 0: get the entry rate
		eq.eM     = p.T/eq.nStar;
		eq.z      = (p.thetaEnt^(1/(1 - p.zetaEnt)))*((1/eq.omegaP)^(p.zetaEnt/(1 - p.zetaEnt)))*((p.zetaEnt)^(p.zetaEnt/(1 - p.zetaEnt)))*...
					((p.delta*eq.VnH(1) + (1 - p.delta)*eq.VnL)^(p.zetaEnt/(1 - p.zetaEnt)));


		%% Step 1: solve firm size distribution and calculate manager demand (xD) and aggregate effort (M)

		[eq,flag] = firmSizeDistBGP(p,eq);

		if flag<0
		    warning('Could not...')
		    diff = ones(3,1)*100000;
		else
		    
			% Manager Demand and Supply
			eq.xD = max((1/p.alpha).*(p.T.*(1/eq.nStar - 1./state)).*eq.muAll,0);
			eq.xD = sum(eq.xD); 
			eq.xS = eq.Ltotal*((p.vartheta/(p.vartheta - 1))*((p.muM*(p.vartheta - 1)/p.vartheta)^p.vartheta)*((eq.omegaM/eq.omegaP)^(p.vartheta-1)));

			% Aggregate effort
			nStarD    = min(alg.NN(end),ceil(eq.nStar));

			partA     = 0;

			for i=1:(nStarD - 1)
			    partA = partA + eq.muAll(i)*(1 - (p.T/i)^p.sigma);
			end
			partB      = (1 - eq.eM^p.sigma)*sum(eq.muAll(nStarD:end));
			eq.M       = ((partA + partB)^(-1));

			% Production worker demand and supply
			eq.lProdD  = 1/(eq.M*eq.omegaP);
			eq.lProdS  = eq.Ltotal*max(0,(1 - ((p.muM*(p.vartheta - 1)/p.vartheta)^p.vartheta)*((eq.omegaM/eq.omegaP)^(p.vartheta))));

			%% Finish
			diff     = [ eq.xD     - eq.xS;...
				   		 eq.lProdD - eq.lProdS;...
			             eq.tauH   - eq.tauHNew];
			           			            
		end 
	end
end

end % function end
