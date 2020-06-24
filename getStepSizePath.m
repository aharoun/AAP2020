function out = getStepSizePath(taut,g,initqUS2qInd,p,alg)

    funcGammaPath = @(t,y)g*y*p.lambda - p.lambda*y*log(y).*interp1(alg.timePath,taut,t);
    solGammaPath  = ode23(funcGammaPath,[alg.timePath2(1) alg.timePath2(end)],(initqUS2qInd^p.lambda)*p.gamma,alg.opt23);
    out           = deval(solGammaPath,alg.timePath);
    out           = [(initqUS2qInd^p.lambda)*p.gamma out];  % to make it consistent in terms of timing.

end
