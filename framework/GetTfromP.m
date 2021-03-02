function [TS_T_ann_interim, extra_c] = GetTfromP(TS_P_ann_interim, HistoricData, info)

    % created 08/10/2019 by Keirnan Fowler, University of Melbourne
    
    % Creates a synthetic annual timeseries of temperature (T) for each 
    % region based on (a) the synthetic timeseries of precipitation P; 
    % (b) the historic relationship between T and P; and (c) a 
    % timeseries of serially correlated random numbers* used to generate 
    % residuals so that the result does not *exactly* follow the line of 
    % best fit from (b)
    
    % * note that only a single timeseries of residuals is generated, since
    % it is assumed that the physical processes leading to the anomolies
    % act across all regions (more or less).  
    
    % first, generate the timeseries of random numbers for use across all regions
    RandomNumbers = rand(info.pars.StochRepLen_yrs, 1); % note, serial correlation comes later
    
    % initialise arrays
    all.T_synth = []; all.resid = []; all.resid_std = []; all.m = []; all.c = []; 

    % for each subarea
    SubareaList = info.SubareaList'; 
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"
        
        % generate a relationship between historic P and T
        P = HistoricData.precip_annual.(subarea{:});
        T = HistoricData.T_annual.(subarea{:});
        [m, c] = LinReg_get_m_and_c(P, T);
        
        % generate a timeseries of model error and calculate stats
        T_mod = m * P + c; 
        resid = T_mod - T; 
        resid_std = std(resid);
        resid_autocorr = kf_autocorr(resid); % disp(['The timeseries of residuals for ' region{1} ' exhibits an autocorrelation of ' num2str(resid_autocorr)])
        
        % generate synthetic rainfall values that perfectly honour the line of best fit
        P_synth = TS_P_ann_interim.(subarea{:});
        T_synth = m * P_synth + c; 
        
        % subject sythetic values to serially correlated random deviations
        Zstd = norminv(RandomNumbers);
        RndDev(1) = 0;
        for iYear = 2:info.pars.StochRepLen_yrs
            AutoCorr = 0.32; % note, a lag-1 autocorrelation of 0.32 is the average value across the 7 regions (see line 31)
            RndDev(iYear) = RndDev(iYear-1)*AutoCorr + resid_std*Zstd(iYear); 
        end
        
        % add modelled values to deviates
        T_synth = T_synth + RndDev'; 
        
        % append variables to arrays
        all.T_synth(:, end+1) = T_synth; all.resid(:, end+1) = resid; all.resid_std(:, end+1) = resid_std; all.m(:, end+1) = m; all.c(:, end+1) = c;         
        
    end
    
    % finalise outputs
    TS_T_ann_interim  = array2table(all.T_synth,   'VariableNames', info.SubareaList); 
    extra_c.resid     = array2table(all.resid,     'VariableNames', info.SubareaList); 
    extra_c.resid_std = array2table(all.resid_std, 'VariableNames', info.SubareaList); 
    extra_c.m         = array2table(all.m,         'VariableNames', info.SubareaList); 
    extra_c.c         = array2table(all.c,         'VariableNames', info.SubareaList); 
    extra_c.RandomNumbers = RandomNumbers; 
end

function autocorr = kf_autocorr(TS)
    temp = corrcoef(TS(1:(end-1)), TS(2:end));
    autocorr = temp(1, 2); 
end

function [m, c] = LinReg_get_m_and_c(x, y)
    P = polyfit(x,y,1);
    m = P(1); c = P(2);
end