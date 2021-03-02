function [TS_PET, extra_g] = GetPETfromT(TS_T, HistoricData, info)
    
    % Created 08/10/2019 by Keirnan Fowler, University of Melbourne
    
    % This code produces, for each region, a monthly timeseries of 
    % synthetic PET based on the monthly timeseries of temperature.  The
    % estimation is based on a regression relationship defined for each 
    % region separately and for each month within each region separately.  
    
    % for reference, see here:
    % ...\01_StochasticClimateGeneration\01_InitialDevelopment\07_PET_from_T\kf_PETfromT_monthly_annotated.png
    
    % get the months of history and also the synthetic data
    months_hist = HistoricData.PET_monthly.Month;
    months_synth = TS_T.Month;
    
    % for each subarea
    SubareaList = info.SubareaList'; TS_PET = [TS_T.Year TS_T.Month]; pars_m = []; pars_c = [];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"
        
        % define vectors with required info
        T_hist = HistoricData.T_monthly.(subarea{:});
        PET_hist  = HistoricData.PET_monthly.(subarea{:});
        T_synth = TS_T.(subarea{:});
        
        % initialise output vector
        PET_synth = nan(size(T_synth));
        
        for iMonth = 1:12
            
            % find rows corresponding to this calendar month
            ind1 = months_hist == iMonth; 
            
            [m, c] = LinReg_get_m_and_c(T_hist(ind1), PET_hist(ind1));
            % scatter(T_hist(ind1), PET_hist(ind1)); hold on; plot([0 40], [c, m*40+c]);
            
            % apply equation to calculate synthetic PET
            ind2 = months_synth == iMonth;
            PET_synth(ind2)  = m.*T_synth(ind2) + c;
            
            % remember the parameters
            m_all(iMonth) = m; 
            c_all(iMonth) = c; 
        
        end       
        
        % append data to arrays
        TS_PET(:, end+1) = PET_synth; 
        pars_m(:, end+1) = m_all; 
        pars_c(:, end+1) = c_all; 
        
    end
    
    % finalise output as tables
    TS_PET = array2table(TS_PET, 'VariableNames', [{'Year'; 'Month'}; info.SubareaList]);
    extra_g.pars_m = array2table(pars_m, 'VariableNames', info.SubareaList); 
    extra_g.pars_c = array2table(pars_c, 'VariableNames', info.SubareaList); 
end

function [m, c] = LinReg_get_m_and_c(x, y)
    P = polyfit(x,y,1);
    m = P(1); c = P(2);
end