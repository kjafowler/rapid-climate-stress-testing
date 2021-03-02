function TS_Q_ML_ForModel = ConvertFlows(TS_Q_subareas_RepCatch_mm, info)
    
    % Created 06/04/2020 by Keirnan Fowler

    % This function takes the flows for the representative catchments of 
    % each subarea and does two things: 
    % - converts from mm/d to ML/d; and
    % - applies pre-specified multiplication factors to calculate the 
    %   final inflows for each point where flow is required.  
    % In this way, a given inflow timeseries can actually be sourced from
    % multiple subareas (eg. see article, Figure #).    
    
    %% for the representative catchment of each subarea, convert from mm/month to ML/month
    RepCatchList = info.FlowConversionFactors.RepCatch'; FlowArray_ReprCatch = [];
    for repcatch = RepCatchList % equivalent of VBA "For Each subarea in SubareaList"
        
        % get flows in mm/month
        ThisFlow_mm = TS_Q_subareas_RepCatch_mm.(repcatch{:});
        
        % get catchment area
        pos = find(strcmp(repcatch, info.RepCatchDetails.identifier));
        CA_km2 = info.RepCatchDetails.CA_km2(pos);
        
        % conversion factor from mm/d to ML/d
        multiplier = 1;                        % mm/d
        multiplier = multiplier/1000;              % m/d
        multiplier = multiplier * CA_km2 * 10^6;   % m³/d or KL/d
        multiplier = multiplier / 1000;            % ML/d
        
        % convert flows to ML/d
        FlowArray_ReprCatch(:, end+1) = ThisFlow_mm*multiplier; 
        
    end
    
    %% for each model input, derive flows as a weighted average of the flows from the different representative catchments
    
    % initialise array to hold output
    TS_Q_ML_ForModel = [TS_Q_subareas_RepCatch_mm.Year TS_Q_subareas_RepCatch_mm.Month];
    
    % derive flows
    ModelInputNames = info.FlowConversionFactors.Properties.VariableNames(2:end); 
    for ModelInput = ModelInputNames % equivalent of VBA "For Each region in RegionList"
        
        % get the factors that apply for this model input
        pos = find(strcmp(ModelInputNames, ModelInput{:})); 
        factors = info.FlowConversionFactors.(ModelInput{:})';
        
        % apply the factors to calculate this timeseries
        ThisFlowArray = FlowArray_ReprCatch.*factors;
        TS_Q_ML_ForModel(:, end+1) = sum(ThisFlowArray, 2);
        
    end
    
    % finalise output as table
    TS_Q_ML_ForModel = array2table(TS_Q_ML_ForModel, 'VariableNames', [{'Year'; 'Month'}; info.FlowConversionFactors.Properties.VariableNames(2:end)']); 
    
end