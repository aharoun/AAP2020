function thetaEnt = getThetaEnt(params,eq)

    if ~isstruct(params)
        p = putParams(params);
    else
        p = params;
    end

    thetaEnt = (eq.z^(1 - p.zetaEnt))*(eq.omegaP^p.zetaEnt)*((p.zetaEnt)^(-p.zetaEnt))*((p.delta*eq.VnH(1) + (1 - p.delta)*eq.VnL)^(-p.zetaEnt));

end
