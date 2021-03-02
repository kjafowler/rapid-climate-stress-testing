function [TS_Q_mm, TS_SM, extra_i] = GetQ_RunRRmodel(TS_P, TS_PET, info)
    
    % Created 04/02/2020 by Keirnan Fowler, University of Melbourne
    
    % Run a pre-calibrated WAPABA monthly rainfall runoff model in order to
    % generate a timeseries of flow (Q, mm/month) for each region. 
    
    % Notes: 
    % The Q timeseries corresponds *not* to the whole-of-region flow but
    % rather to flow from a selected single catchment from that region.
    % In each case the catchment has been selected because it is the 
    % one with the most has similar climatology to the region as a whole,
    % compared to other available HRS catchments in the region (HRS stands 
    % for 'Hydrologic Reference Station' and is the BoM's list of high
    % quality flow gauging stations).
    
    % In each of the 7 cases, the selected catchment has previously been 
    % calibrated using region-wide climatic data (specifically, monthly 
    % timeseries of region-average Precip. and region-average PET).  The
    % calibration method was multi-objective optimisation using AMALGAM,
    % where one objective is performance over a dry period (7 driest years)
    % and the other objective is performance over the remainder of the
    % timeseries.  This method is as per Keirnan's PhD (see here
    % https://doi.org/10.1002/2015WR018068).  The selected parameter set
    % was the one that was closest to the 'perfect point' in this
    % 2-dimensional space.  
    
    % Calibrated in this way, the models attained seemingly good scores -
    % most cases exceeded KGE of 0.9 and almost all exceeded KGE of 0.8
    % except the exceptionally dry Region G for which [KGEnondry, KGEdry] =
    % [0.85, 0.6].  Furthermore, the degradation in performance due to using 
    % regional climatic data was minimal, see 
    % /data/cephfs/punim0002/kfowler/analysis/2020/20200203_GoulburnLP_WAPABAcal_RegionalForcingData/outputs/ResultsComparison_FinishingTouches.png
    
    % On this basis, it is assumed the models are fit to be applied in the
    % decision scaling runs driven only by regional-scale climatic inputs.
    % It is noted that it is *not* assumed that the rainfall-runoff
    % relationship (as described by WAPABA) is stationary: this is
    % perturbed separately at a later stage in the code.  
    % 
    % Note, we later derive scaling factors to relate the catchment-scale 
    % model outputs to the river model inflow points (not done here).  
    
    % for each subarea
    SubareaList = info.SubareaList'; TS_Q_mm = [TS_P.Year TS_P.Month]; TS_SM = [TS_P.Year TS_P.Month];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"
        
        % define inputs
        P            = TS_P.(subarea{:});
        PET          = TS_PET.(subarea{:});
        ListOfMonths = TS_P.Month;
        
        % get WAPABA parameters from prior calibration
        WAPABApars = info.WapabaParSets.(subarea{:});
        
        % run WAPABA
        ParInput.a1   = WAPABApars(1, 1); 
        ParInput.a2   = WAPABApars(1, 2); 
        ParInput.B    = WAPABApars(1, 3); 
        ParInput.Smax = WAPABApars(1, 4); 
        ParInput.K    = WAPABApars(1, 5); 
        
        [SimFlows, SimVars] = RunWapaba(ParInput, P, PET, ListOfMonths);
        
        % provide a more appropriate name that reflects both the region and the catchment name
        pos = find(strcmp(subarea, info.SubareaList));
        Qname = info.RepCatchDetails.identifier(pos);
        
        % append data to arrays
        TS_Q_mm(:, end+1) = SimFlows; TS_SM(:, end+1) = SimVars.States.S'/ParInput.Smax;
        extra_i.(subarea{:}).SimVars = SimVars; 
    end
    
    % finalise output as tables (note use of different identifiers for representative catchments) 
    TS_Q_mm = array2table(TS_Q_mm, 'VariableNames', [{'Year'; 'Month'}; info.RepCatchDetails.identifier]);
    TS_SM   = array2table(TS_SM,   'VariableNames', [{'Year'; 'Month'}; info.RepCatchDetails.identifier]);
    
end