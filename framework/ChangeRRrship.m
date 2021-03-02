function TS_Q_shifted = ChangeRRrship(TS_Q_unpert, deltaRRrship, info)
    
    % apply a shift to the annual rainfall runoff relationship.  
    
    % after Saft et al., 2015, we apply this shift in Box-Cox transformed
    % space.  This code assumes that a single Box-Cox parameter is applicable
    % across all subareas.  Supp. Mat. ## shows that this is a valid 
    % assumption for the case study used in the paper but this does not 
    % guarantee it is a valid assumption elsewhere, so please consider this 
    % carefully.  
    
    % In any case, the perturbation in RR relationship is defined as a 
    % vertical shift in the transformed space.  Thus, the shifts can only
    % be understood with reference to the transformed space they are defined
    % in.  
    
    % set month to start on, based on definition of water years
    month_vector = TS_Q_unpert.Month; 
    months = {'January', 'February', 'March', 'April', 'May', 'June', 'July', ...
          'August', 'September', 'October', 'November', 'December'};
    wy_StartMonth = find(strcmp(months, info.pars.WaterYearStart_flow)); 
    
    % find first complete year
    cal_year_vector = TS_Q_unpert.Year; 
    pos = find(month_vector == wy_StartMonth, 1, 'first'); % first instance of starting month
    FirstYear = cal_year_vector(pos); first_pos = pos; 
    
    % find last complete year
    pos = find(month_vector == wy_StartMonth, 1, 'last'); % last instance of starting month
    NumMonthsLeft = size(TS_Q_unpert, 1) - pos + 1; 
    if NumMonthsLeft < 12, pos = pos - 12; end
    LastYear = cal_year_vector(pos); 
    NumWY = LastYear - FirstYear + 1; 
    
    % work with array version of the data (faster)
    TS_Q_array = table2array(TS_Q_unpert); 
    TS_Q_array = TS_Q_array(:, 3:end); % delete cols for Year and Month
    
    % for each subarea
    SubareaList = info.RepCatchDetails.identifier'; TS_Q_shifted_array = nan(size(TS_Q_array)); iCol = 0; 
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"    
        
        iCol = iCol + 1; 
        
        pos_endmonth = first_pos - 1; 
        for iYear = 1:NumWY

            % get index: which rows correspond to months of the current year?
            pos_startmonth = pos_endmonth + 1; pos_endmonth = pos_startmonth + 11; 
            ind = pos_startmonth:pos_endmonth;
            
            % calculate the annual flow
            ThisYear_MonthlyQ = TS_Q_array(ind, iCol);
            ThisYear_Q = sum(ThisYear_MonthlyQ); 
            
            % Box Cox transform
            lambda = info.pars.Streamflow_BoxCox_Lambda; % see Section ## / Supp Mat ##. 
            Q_BC = (ThisYear_Q^lambda - 1)/lambda; 
            
            % shift in rr relationship
            Q_BC_shifted = max(0, Q_BC + deltaRRrship); % see notes above
            
            % reverse Box Cox transform
            Q_shifted = ((Q_BC_shifted * lambda)+1)^(1/lambda);
            
            % disaggregate to monthly flow using same pattern as before
            weights = ThisYear_MonthlyQ / ThisYear_Q; 
            Shifted_MonthlyQ = weights*Q_shifted; 
            
            % place back in array
            TS_Q_shifted_array(ind, iCol) = Shifted_MonthlyQ; 
        end
    end
    
    % transfer data from array to output table
    TS_Q_shifted = array2table([cal_year_vector month_vector TS_Q_shifted_array], 'VariableNames', [{'Year'; 'Month'}; info.RepCatchDetails.identifier]);
    
end