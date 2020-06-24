
function [lambdaModel,qUS2qINPathModel,timePath] = calibrateDiffusion(pIN,eqIN,pUS,eqUS)


pIN.gamma = pUS.gamma;

fid = fopen('Moments/ctfpINpwt9.csv');
C   = textscan(fid, '%f%f','delimiter',',');
fclose(fid);

yearEnd  = 2005;
yearInit = 1985;
tfpUS2tfpINPathData = 1./C{2}((C{1} <= yearEnd)&(C{1} >= yearInit));
qUS2qINPathData     = (eqIN.lProdS/eqUS.lProdS)*(eqIN.M/eqUS.M)*tfpUS2tfpINPathData;
timePath            = linspace(0,yearEnd-yearInit,yearEnd-yearInit+1);           % for calibration

report = 0;
options    = optimset('Display','off');
lambdaModel = fminsearch(@(x)distanceZ(x),[.1,2],options);

timePath              = linspace(0,2400-yearInit,2400-yearInit+1);              % for reporting
report = 1;
[~,qUS2qINPathModel] = distanceZ(lambdaModel);                                  % calibrate both lambda and z0

function [dist,qUS2qINPathModel] =  distanceZ(x)

    in = x(1);
    z0 = x(2);

    a = eqUS.g - log(pIN.gamma)*eqIN.tau ;
    b = in*eqIN.tau;
    c = z0;

    qUS2qINPathModel  = exp((a - exp(-b*timePath)*(a - b*log(c)))./b);

    if report == 0
        dist = sum((qUS2qINPathData - qUS2qINPathModel').^2);
    else
        dist = inf;
    end


end


end


