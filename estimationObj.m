function [score,pIN,eqIN,mIN,pUS,eqUS,mUS,meanReg] = estimationObj(parin,vers)
    % Objective function

    global alg scoreBest

    nan_score  = Inf;
    inf_score  = Inf;
    bnd_score  = Inf;
    fail_score = Inf;


    ptagIN = ['India' vers];
    ptagUS = ['US' vers];
    
    initAlg();
    [pIN,pNamesIN] = readParam(['params' filesep 'params-' ptagIN '.txt']);
    [pUS,pNamesUS] = readParam(['params' filesep 'params-' ptagUS '.txt']);

    
    if any(strfind(ptagIN,'Payment')) 
        estParams   = [1:7];                % estimate sigma
        bloomMoment = false;                % don't use Bloom moment 
        if any(strfind(ptagIN,'BloomAndPayment'))
            bloomMoment = true;             % use Bloom moment
        end
    else 
        estParams   = [1:2 4:7];            % do not estimate sigma
        pIN.sigma = parin(end)*pIN.sigma;   % in this case, last element of parin is for sigma
        pUS.sigma = pIN.sigma;              % common sigma across countries
        bloomMoment = true;                 % use Bloom moment
    end    
        
        
    iter = 1;
    for i = 1:length(estParams)
        pIN.(pNamesIN{estParams(i)}) = parin(iter)*pIN.(pNamesIN{estParams(i)});
        iter = iter + 1;
    end

    for i = 1:length(estParams)
    	pUS.(pNamesUS{estParams(i)}) = parin(iter)*pUS.(pNamesUS{estParams(i)});
    	iter = iter + 1;
    end

    % Identification of muM
    %-----------------------------------------------------------------------------
    initAlg(ptagUS); 
    fid = fopen(alg.momentDataFile);
    C   = textscan(fid, '%s%f','delimiter',',');
    fclose(fid);
    dataUS = cell2struct(num2cell(C{2}),C{1});

    initAlg(ptagIN); 
    fid = fopen(alg.momentDataFile);
    C   = textscan(fid, '%s%f','delimiter',',');
    fclose(fid);
    dataIN = cell2struct(num2cell(C{2}),C{1});
    
    % With this, we hit relative managerial share of Indian migrants exactly
    pIN.muM   = ((dataIN.shareMan/dataIN.shareManMigrantIN)*(dataIN.shareManMigrantUS/dataUS.shareMan))^(1/pUS.vartheta);
    %-----------------------------------------------------------------------------

    % parameter limits
    paramsIN = (cell2mat(struct2cell(pIN)))';
    paramsUS = (cell2mat(struct2cell(pUS)))';

    plbIN            = zeros(1,length(paramsIN)) + 0.00000000001;    
    pubIN            = Inf*ones(1,length(paramsIN));  
    pubIN(3)         = .99999999999;
    plbIN(11)		 = -10; 	

    plbUS            = zeros(1,length(paramsUS)) + 0.00000000001;    
    pubUS            = Inf*ones(1,length(paramsUS));  
    pubUS(3)         = .99;
    plbUS(11)		 = -10; 	


    if any(paramsIN < plbIN) || any(paramsIN > pubIN) || any(paramsUS < plbUS) || any(paramsUS > pubUS) 
    	fprintf(1,'BOUNDS ERROR\n');
    	score = bnd_score;
    else
    	% INDIA
    	alg = struct;
    	initAlg(ptagIN);
    	[eqIN,pIN,flagIN] = solveBGP(pIN, alg.cnt); 
    	if (flagIN == 0)
    		fprintf(1,'EQ SOLVE FAILED\n');
    		score = fail_score;
    	else
    		[~,mIN,eqIN] = callMomentBGP(pIN,eqIN,0);
    		
    		% US
    		alg = struct;
    		initAlg(ptagUS);
    		[eqUS,pUS,flagUS] = solveBGP(pUS, alg.cnt);
    		if (flagUS == 0)
    			fprintf(1,'EQ SOLVE FAILED\n');
    			score = fail_score;
    		else
    			[~,mUS,eqUS] = callMomentBGP(pUS,eqUS,0);

    			% make sure step size is consistent with 0.02 growth
    			eqUS.g = 0.02;
    			pUS.gamma = eqUS.stepSizeSS;
    			pIN.gamma = pUS.gamma;
    			
    			% calculate catch-up parameters
    			[lambdaModel,~,~] = calibrateDiffusion(pIN,eqIN,pUS,eqUS);
    			pIN.lambda = lambdaModel(1);
                
                if bloomMoment 
                    %-------------------------------------------------------------------------------------
                    % Bloom Moment
                    pool = gcp('nocreate');
                    if isempty(pool)
                        parpool('local',alg.nPool);
                    end
    			    bloomWeight = 1.0;
                    [meanReg,~,~,flagBloom] = bloomExercise(eqIN,pIN,eqUS);
                else
                    flagBloom   = 1;
                    bloomWeight = 0.0;
                    meanReg     = 1000.0;
                end

                %-------------------------------------------------------------------------------------  
                    % Now put things together
                momentWgt   = [mIN.momentWgt    mUS.momentWgt	bloomWeight];  % the last one is for Bloom moment
            	mErr        = [mIN.mErr	        mUS.mErr        dataIN.bloomMoment - meanReg];
            	momentModel = [mIN.momentModel  mUS.momentModel  meanReg];
            	momentData  = [mIN.momentData   mUS.momentData dataIN.bloomMoment];
                
            	score = (sum((momentWgt'.*(abs(mErr')./(0.5*abs(momentModel') + 0.5*abs(momentData'))))))/length(momentModel);

                if flagBloom==0 
                    fprintf(1,'EQ SOLVE FAILED AT BLOOM EXERCISE\n');
                    score = fail_score;
                end
    								
        		if score<=scoreBest
        	   		initAlg(ptagIN);
        			writeOutput(pIN,mIN,'logs/outputBestIN','');
        			initAlg(ptagUS);
        			writeOutput(pUS,mUS,'logs/outputBestUS','');
        			scoreBest = score;
        			fprintf('BLOOM TARGET = %2.3f, BLOOM MODEL = %2.8f, sigma = %2.8f\n',dataIN.bloomMoment, meanReg, pIN.sigma)
        			fprintf(1,'Score =\t %f\n',score);
    		    end
        	end
    	end
    end
end


