function [eq,m,p]= solver(cnt)
    % Solves the equilibrium of the model and calculate the moments.
    global alg
    alg = struct;
    if nargin == 0
        initAlg();
        fprintf('Country: %s\n',alg.ssub)
    else
    	initAlg(cnt);
        fprintf('Country: %s\n',cnt)        
    end 

    [p, ~]     = readParam(alg.paramsFile);    
    [eq, p, ~] = solveBGP(p,cnt);
    [~, m, eq] = callMomentBGP(p,eq,1);
end
