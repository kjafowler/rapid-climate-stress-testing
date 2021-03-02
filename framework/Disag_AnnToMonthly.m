function [TS_P_interim, TS_T, extra_f] = Disag_AnnToMonthly(TS_P_ann, TS_T_ann, HistoricData, info)
    
    % created 08/10/2019 by Keirnan Fowler, University of Melbourne
    
    % Disaggregate P and T timeseries to monthly by method of fragments
    % Do this by randomly choosing a year from history and copying:
    % (i)  its proportions of rainfall in each year
    % (ii) its seasonal pattern of temperature (but shift it laterally, 
    %      since we've already perturbed the annual temp values).
    
    % reset the random seed to ensure the same sequence of random number is used at each position in the space
    if ~info.isOctave, rng default; end
    
    % choose a year from history for each synthetic year
    option = 3; % uses wet, medium and dry 
    years = RandomlyAssignYears(HistoricData, TS_P_ann, info, option);
    
    % for each subarea
    SubareaList = info.SubareaList'; TS_P_interim = []; TS_T = [];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"
            
        % intialise input vectors (makes code faster)
        monthly_P_hist_all    = HistoricData.precip_monthly.(subarea{:}); 
        monthly_T_hist_all    = HistoricData.T_monthly.(subarea{:}); 
        monthly_P_ListOfYears = HistoricData.precip_monthly.Year; 
        this_TS_P_ann         = TS_P_ann.(subarea{:}); 
        this_TS_T_ann         = TS_T_ann.(subarea{:}); 

        % intialise output vectors (makes code faster)
        P_synth_mon_all    = nan(1, info.pars.StochRepLen_yrs*12); 
        T_synth_mon_all    = nan(1, info.pars.StochRepLen_yrs*12); 
        year_synth_all      = nan(1, info.pars.StochRepLen_yrs*12); 
        month_synth_all     = nan(1, info.pars.StochRepLen_yrs*12); 

        % loop over years
        iMonth = 1; 
        for iYear = 1:info.pars.StochRepLen_yrs

            ind = monthly_P_ListOfYears == years(iYear);

            % precip
            monthly_P_hist = monthly_P_hist_all(ind);
            pattern = monthly_P_hist / sum(monthly_P_hist);
            P_synth_ann = this_TS_P_ann(iYear);
            P_synth_mon = P_synth_ann*pattern;
            
            % temperature
            monthly_T_hist = monthly_T_hist_all(ind);
            averageT_hist = sum(monthly_T_hist .* [31 28.25 31 30 31 30 31 31 30 31 30 31]')/365.25;
            T_synth_ann = this_TS_T_ann(iYear);
            T_synth_mon = monthly_T_hist + (T_synth_ann-averageT_hist);

            % add to vectors
            P_synth_mon_all(iMonth:iMonth+12-1) = P_synth_mon(1:12); 
            T_synth_mon_all(iMonth:iMonth+12-1) = T_synth_mon(1:12); 
            year_synth_all  (iMonth:iMonth+12-1) = repmat(iYear, 1, 12); 
            month_synth_all (iMonth:iMonth+12-1) = [1 2 3 4 5 6 7 8 9 10 11 12];
            % subplot(1, 2, 1); plot(P_synth_mon_all(1:iMonth+12-1), 'b'); 
            % subplot(1, 2, 2); plot(T_synth_mon_all(1:iMonth+12-1), 'r'); 
            iMonth = iMonth + 12; 

        end
        
        % append to array
        TS_P_interim(:, end+1) = P_synth_mon_all;        
        TS_T        (:, end+1) = T_synth_mon_all; 
    end
    
    % include year and month in arrays
    TS_P_interim = [year_synth_all' month_synth_all' TS_P_interim];
    TS_T         = [year_synth_all' month_synth_all' TS_T];
    
    % finalise output as tables
    TS_P_interim = array2table(TS_P_interim, 'VariableNames', [{'Year'; 'Month'}; info.SubareaList]); 
    TS_T         = array2table(TS_T,         'VariableNames', [{'Year'; 'Month'}; info.SubareaList]); 
    
    % add ancillary information to data structure extra_f
    extra_f.years = years; 
    
end

function years = RandomlyAssignYears(HistoricData, TS_P_ann, info, option)
    
    lastyear = max(HistoricData.precip_monthly.Year);
    firstyear = min(HistoricData.precip_monthly.Year);
    
    % define 'indicator series': a representative series that is used to
    % categorise year by year into 'wet', 'medium' or 'dry'
    SubareaList = info.SubareaList;     
    IndSer_hist  = HistoricData.precip_annual.(SubareaList{1}); % assume subarea 1 is sufficient indicator of wet and dry years
    IndSer_synth = TS_P_ann.(SubareaList{1}); 
    
    switch option
        case 1 % option 1: randomly choose a year from history, without 
               % matching wet years to wet years and dry years to dry years
            years = randi([firstyear lastyear],1, info.pars.StochRepLen_yrs);
            
        case 2 % option 2: randomly choose a year from history but in three
               % bins: dry, medium, and wet. Divide synthetic years
               % into these bins so that an equal number go in each bin.
               % Then do the same with history. Each synthetic year is only
               % matched with historic years from its bin.  
            [div3_hist,  ind_h] = DivideIntoThree(IndSer_hist);
            [div3_synth, ind_s] = DivideIntoThree(IndSer_synth);
            for iYear = 1:info.pars.StochRepLen_yrs
                switch div3_synth{iYear}
                    case 'dry',    AvailableYears = find(ind_h == 1); 
                    case 'medium', AvailableYears = find(ind_h == 2); 
                    case 'wet',    AvailableYears = find(ind_h == 3); 
                end
                pos = randi(size(AvailableYears, 2)); years(iYear) = AvailableYears(pos); 
            end
            years = years + firstyear - 1; 
            
        case 3 % option 3: same as option 2 except that the thresholds from 
               % history are considered to apply to syntethic future data
               % also.  This means that (eg) if we are investigating a future
               % reduction in rainfall, more than 1/3 of synthetic years
               % will be drawn from the 'dry' bin (and vice versa).  This
               % implicitly means that the seasonality will change in all
               % cases where the deltaP axis value is not zero.  
            [div3_hist, ind_h, threshold_hist] = DivideIntoThree(IndSer_hist);
            for iYear = 1:info.pars.StochRepLen_yrs
                if IndSer_synth(iYear) <= threshold_hist(1)
                    i = 1; ind_s(iYear) = i; AvailableYears = find(ind_h == i); 
                elseif IndSer_synth(iYear) <= threshold_hist(2)
                    i = 2; ind_s(iYear) = i; AvailableYears = find(ind_h == i); 
                else
                    i = 3; ind_s(iYear) = i; AvailableYears = find(ind_h == i); 
                end
                pos = randi(size(AvailableYears, 2)); years(iYear) = AvailableYears(pos); 
            end
            years = years + firstyear - 1; 
    end
end

function [div3_hist, ind, threshold_hist] = DivideIntoThree(IndSer)

    % set thresholds
    threshold_hist(1) = prctile(IndSer, 33.333); % upper bound dry, lower bound medium
    threshold_hist(2) = prctile(IndSer, 66.666); % upper bound medium, lower bound wet
    
    % categorise each year
    for iYear = 1:size(IndSer, 1)
        if IndSer(iYear) <= threshold_hist(1)
            div3_hist{iYear} = 'dry';    ind(iYear) = 1; 
        elseif IndSer(iYear) <= threshold_hist(2)
            div3_hist{iYear} = 'medium'; ind(iYear) = 2; 
        else
            div3_hist{iYear} = 'wet';    ind(iYear) = 3; 
        end
    end

end