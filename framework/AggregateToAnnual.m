function HistoricData = AggregateToAnnual(HistoricData, pars)

    % created January 2021 by Keirnan Fowler, University of Melbourne
    
    % Format data by aggregating monthly to annual.  This function was
    % originally designed for aggregating historic data (that's why it has
    % checks for missing data and also explains the variable names) but is 
    % also used for aggregating stochastic data.  
    
    for iDataType = 1:5 % P, T, PET, and two types of Q - see lists below
        
        % initialise variables
        [data_monthly, ... % extracted data just for this data type
         agg_type, ...     % whether we are summing over the year, or taking the mean
         wy, ...           % vector of water years (not yet vetted for incompleteness)
         wy_list] ...      % water year for each row of data_monthly
         = init_agg(HistoricData, pars, iDataType); 
        
        % check which years have all 12 values and find the earliest and latest years to include
        [wy_complete, year_range] = CheckData(data_monthly, wy, wy_list, pars); 

        % undertake aggregation to annual, writing the results to an array
        data_annual_array = aggregate(data_monthly, agg_type, wy, wy_list, wy_complete); 
        
        % format as table and clip to active years
        data_annual = ToTable(data_annual_array, data_monthly, wy, year_range);
        
        % pass to output structure
        switch iDataType
            case 1, HistoricData.precip_annual                            = data_annual; 
            case 2, HistoricData.T_annual                                 = data_annual; 
            case 3, HistoricData.PET_annual                               = data_annual; 
            case 4, HistoricData.flow_annual.ReachInflows_ML              = data_annual;              
            case 5, HistoricData.flow_annual.RepresentativeCatchments_mm  = data_annual; 
        end        
        
        clearvars -except HistoricData pars iDataType
        
    end
end

%% functions

function [data_monthly, agg_type, wy, wy_list] = init_agg(HistoricData, pars, iDataType)
        
    switch iDataType
        case 1, data_monthly = HistoricData.precip_monthly;                           agg_type = 'sum' ; wy_start_txt = pars.WaterYearStart_clim; 
        case 2, data_monthly = HistoricData.T_monthly;                                agg_type = 'mean'; wy_start_txt = pars.WaterYearStart_clim; 
        case 3, data_monthly = HistoricData.PET_monthly;                              agg_type = 'sum' ; wy_start_txt = pars.WaterYearStart_clim; 
        case 4, data_monthly = HistoricData.flow_monthly.ReachInflows_ML;             agg_type = 'sum' ; wy_start_txt = pars.WaterYearStart_flow; 
        case 5, data_monthly = HistoricData.flow_monthly.RepresentativeCatchments_mm; agg_type = 'sum' ; wy_start_txt = pars.WaterYearStart_flow; 
    end
    cal_year_vector = data_monthly.Year; month_vector = data_monthly.Month;
    
    % set month to start on, based on definition of water years
    months = {'January', 'February', 'March', 'April', 'May', 'June', 'July', ...
          'August', 'September', 'October', 'November', 'December'};
    wy_StartMonth = find(strcmp(months, wy_start_txt)); 
    
    % find first complete year
    pos = find(month_vector == wy_StartMonth, 1, 'first'); % first instance of starting month
    FirstYear = cal_year_vector(pos); first_pos = pos; 
    
    % find last complete year
    pos = find(month_vector == wy_StartMonth, 1, 'last'); % last instance of starting month
    NumMonthsLeft = size(data_monthly, 1) - pos + 1; 
    if NumMonthsLeft < 12, pos = pos - 12; end
    LastYear = cal_year_vector(pos); 
    NumWY = LastYear - FirstYear + 1; 
    wy = FirstYear:LastYear; 
    
    % create wy_list, a list of water years for each month in the historic record
    wy_list = nan(1, size(data_monthly, 1)); 
    pos_endmonth = first_pos - 1; 
    for iYear = 1:NumWY
        
        % get index: which rows correspond to months of the current year?
        pos_startmonth = pos_endmonth + 1; pos_endmonth = pos_startmonth + 11; 
        ind = pos_startmonth:pos_endmonth;
        
        % for each month in this water year, set wy_list to the calendar
        % year of the first month in the water year... got that? ;-) 
        wy_list(ind) = cal_year_vector(pos_startmonth); 
        
    end
end

function [wy_complete, year_range] = CheckData(data_monthly, wy, wy_list, pars)
        
    % initialise variables
    NumCols = size(data_monthly, 2) - 2; % don't include year, month
    NumYears = size(wy, 2); 
    
    % add default preferences in the case where pars.checking is not specified
    if ~isfield(pars, 'checking')
        if NumYears < 200 % assume historic data
            pars.checking = true; 
        else % assume stochastic data - no gaps
            pars.checking = false;
        end
    end
    
    if ~pars.checking % no need to check
        wy_complete(1:NumYears, 1:NumCols) = true;
        year_range.max = max(wy);
        year_range.min = min(wy);
    else
        
        % more initialisation
        wy_complete(1:NumYears, 1:NumCols) = false; % initially assume incomplete years
        year_range.max = min(wy); % initialise as arbitrarily early year
        year_range.min = max(wy); % initialise as arbitrarily late year
        data_monthly_array = table2array(data_monthly);

        % for each column, check which years have all 12 values
        for iCol = 1:NumCols
            for iYear = 1:NumYears
                ThisYear = wy(iYear);
                ind = wy_list==ThisYear; 
                vals = data_monthly_array(ind, iCol+2); 
                if nnz(~isnan(vals)) == 12
                    wy_complete(iYear, iCol) = true; 
                    if ThisYear > year_range.max, year_range.max = ThisYear; end
                    if ThisYear < year_range.min, year_range.min = ThisYear; end
                end
            end
        end
    end
end

function data_annual_array = aggregate(data_monthly, agg_type, wy, wy_list, wy_complete)
    NumCols = size(data_monthly, 2) - 2; % don't include year, month
    data_monthly_array = table2array(data_monthly);
    NumYears = size(wy, 2); 
    data_annual_array(1:NumYears, 1:NumCols) = nan;
    for iCol = 1:NumCols
        for iYear = 1:NumYears
            ThisYear = wy(iYear); 
            if iYear == 1, ind = find(wy_list==ThisYear); else, ind = (ind(12)+1):ind(12)+12; end
            vals = data_monthly_array(ind, iCol+2); 
            if wy_complete(iYear, iCol)
                switch agg_type
                    case 'sum',  data_annual_array(iYear, iCol) = sum(vals); 
                    case 'mean', data_annual_array(iYear, iCol) = mean(vals); 
                end
            end
        end
    end
    %%
end

function data_annual = ToTable(data_annual_array, data_monthly, wy, year_range)
    
    % get data and clip as appropriate
    sr = find(wy == year_range.min); % start row
    er = find(wy == year_range.max); % end row
    array_for_table = [wy(sr:er)' data_annual_array(sr:er, :)];
    
    % get headers
    orig_var_names = data_monthly.Properties.VariableNames;
    new_var_names = orig_var_names(1, [1 3:end]); % omit column 2 (month)
    
    % create table
    data_annual = array2table(array_for_table, 'VariableNames', new_var_names); 
    
end