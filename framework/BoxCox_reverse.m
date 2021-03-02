function TS_P_ann_HighFreq = BoxCox_reverse(Matalas_denormalised, BoxCoxParams, info)

    % created 07/10/2019 by Keirnan Fowler, University of Melbourne
    
    % Reverse the earlier box cox transformations
    
    subareas = info.SubareaList'; out_all = [];
    for subarea = subareas  % equivalent of VBA "For Each subarea in subareas"
        
        % get data for this site
        ThisData   = table2array(Matalas_denormalised(:, subarea));
        pos = strcmp(BoxCoxParams.subarea, subarea);
        ThisLambda = BoxCoxParams.BoxCoxParam1(pos);
        ThisShift  = BoxCoxParams.BoxCoxParam2(pos);
        
        % for each element in this vector, reverse the shift
        for iYear = 1:info.pars.StochRepLen_yrs
            ThisData(iYear) = (((ThisLambda * ThisData(iYear))+1)^(1/ThisLambda)) - ThisShift;
        end
        
        out_all = [out_all ThisData]; % accumulate in array
        
    end 
    
    % format into table
    TS_P_ann_HighFreq = array2table(out_all, 'VariableNames', info.SubareaList); 
end