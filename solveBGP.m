function [eq,p,flag1]  = solveBGP(params,cnt)
    % propose omegaM omegaP and tauH and check the labor market conditions and new tau
    global alg

    if ~isstruct(params)
        p = putParams(params);
    else
        p = params;
    end


    init    = load(alg.eqvarsFile);

    if init(1) < p.sigma*p.alpha % initial nStar would be too low!
        init(1) = (p.sigma*p.alpha)*4.0;
    end

    objFnc = @(x)solveBGPFunc(x,p);
    [res,fval,flag] = fsolve(objFnc,init, alg.mopts);

    if flag<=0 || sum(fval.^2)>1.0e-08
        init(3) = init(3)*5.0;
        init(1) = init(1)*2.0;
        func2 = @(x)sum(objFnc(x).^2);
        [res,~,~]    = fminsearch(func2,init,alg.mopts);
        [res,fval,~] = fsolve(objFnc,res,alg.mopts);
    end

    [~,eq] = solveBGPFunc(res,p);

    if  sum(fval.^2)<1.0e-08
            flag1 = 1;
            if alg.saveInit == 1 
                save(alg.eqvarsFile,'res','-ascii','-double');
            end    

            %% additional objects
            if strcmp(cnt(1:2), 'US') 
                eq.g          = log(p.gamma)*eq.tau;          % growth rate
                eq.stepSizeSS = exp(0.02/eq.tau);             % recover step size consistent with 0.02 growth
            elseif strcmp(cnt(1:2), 'In') 
                eq.g          = p.gUS;                        % at BGP, it grows same rate with the US.
                eq.stepSizeSS = exp(eq.g/eq.tau);             % recover the step size at steady state
            end
            
            eq.r_gY = p.rho;
            eq.r    = eq.g + p.rho;                           % interest rate
            eq      = ageSizeDistBGP(p,eq,alg.maxAgeEst + 1); % get age - size distribution

            eq = orderfields(eq);
    else
            disp('Could not solve it!')
            eq    = 0;
            flag1 = 0;
    end

end
