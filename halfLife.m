function out = halfLife(tauIN,gUS,lambda,gamma,z0)

    % half life of qIN2QUS, NOT qUS2QIN    
    a = gUS- log(gamma)*tauIN ;
    b = lambda*tauIN;
    c = z0;
    out = (-1/b)*log((a-b*log(1/(0.5*(1/exp(a/b)+1/c))))/(a-b*log(c)));

end