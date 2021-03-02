function [Precip_BoxCox, BoxCoxParams] = BoxCox(HistoricData, UseBoxCox, info)

    % created 02/10/2019 by Keirnan Fowler, University of Melbourne
    
    % The purpose of this code is to conduct a shifted Box-Cox
    % transformation on the high-frequency component of precipitation.  
    % This is done separately for each subarea.   
    
    % First, some notes about the Box-Cox transformation.  The most common
    % form has a single parameter.  In general, this form would be suitable
    % if we were applying it to unaltered annual precipitation timeseries.
    % However, it cannot be applied to data with negative values.  Here, we
    % are applying it to the high frequency component of annual
    % precipitation, a timeseries which has a mean of zero and thus many
    % negative values.  
    
    % An alternative is the two parameter Box-Cox transformation and this
    % is applied here.  For the two parameter version, a shift is added to 
    % the data to remove negatives.  One drawback of this method is the need 
    % to choose the value of the second parameter.  Here, the data 
    % originated as high-frequency oscillations around a non-zero mean.
    % Given this, it seems logical to set the second 'shift' parameter
    % value to the same value as the mean in the original data (ie the long
    % term average rainfall). This is done individually for each subregion.
    
    if ~UseBoxCox
        % if no Box Cox'ing, just set the first parameter to 1 and the 
        % second parameter to zero (equivalent to no transformation)
        BoxCoxParam1(1:info.pars.NumSubareas) = 1; BoxCoxParam2(info.pars.NumSubareas) = 0; 
        Precip_BoxCox = HistoricData.precip_ann_HighFreq; 
    else
        
        % initialise
        Precip_BoxCox_array = HistoricData.precip_annual.Year;
        SubareaList = info.SubareaList';
        BoxCoxParam1 = []; BoxCoxParam2 = []; 
        
        for subarea = SubareaList % equivalent of VBA "For Each Subarea in Subareas"
            
            ThisData = HistoricData.precip_ann_HighFreq.(subarea{:});
            
            % set the parameter 2 (shift) value to the long term mean (see notes above) 
            shift = mean(HistoricData.precip_annual.(subarea{:}));

            % find the best parameter 1 (lambda) value
            fun = @(lambda)abs(skewness((((ThisData+shift).^abs(lambda))-1)/abs(lambda)));
            lambda0 = 1.0;
            lambda = abs(fminsearch(fun, lambda0)); % choose value that minimises the skewness
            TransformedData = (((ThisData+shift).^abs(lambda))-1)/abs(lambda); 
            
            % store the transformed data and adopted parameters
            Precip_BoxCox_array(:, (end+1)) = TransformedData';
            BoxCoxParam1(end+1) = lambda; BoxCoxParam2(end+1) = shift; 

        end
    end
    
    % create table for output 1
    Precip_BoxCox = array2table(Precip_BoxCox_array, 'VariableNames', ['Year'; info.SubareaList]);
    
    % create table for output 2
    subarea = SubareaList'; BoxCoxParam1 = BoxCoxParam1'; BoxCoxParam2 = BoxCoxParam2'; 
    BoxCoxParams = table(subarea, BoxCoxParam1, BoxCoxParam2); 
end