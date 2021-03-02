function [TS_Q, TS_P, TS_T, TS_PET, TS_SoilMoisture, extra] = GetSingleReplicate_clim(deltaP, deltaT, deltaLowFreqP, deltaSeasonality, HistoricData, info)

    % Created October 2019, updated January 2021 by Keirnan Fowler, University of Melbourne
    
    % Integrated framework for rapid climate stress testing on a monthly timestep
    % by Keirnan Fowler, Natasha Ballis, Avril Horne, Andrew John, Rory Nathan and Murray Peel 

    % Licence: CC BY 3.0 - see https://creativecommons.org/licenses/by/3.0/au/    
    
    % The purpose of this code is generate a single replicate of arbitrary
    % length (typically >10^3 years) of rainfall and temperature/PET data,
    % for a set of multiple sub-areas, replicating key statistics of the 
    % historic data, except perturbed to match the specified perturbations
    % (deltaP, deltaT, deltaLowFreqP and deltaSeasonality).  
    
    % A further output is the timeseries of flow, but note that this flow
    % timeseries has not yet had any perturbation in rainfall-runoff 
    % relationship applied.
    
    % This code consists of the following steps:
    % (a) Generate timeseries (TS) for high frequency component of annual 
    %     precipitation, unperturbed
    % (b) Generate low frequency component of annual precipitation including
    %     perturbation of Hurst coefficient (low frequency behaviour), and
    %     add the high and low frequency components together
    % (c) Generate timeseries of annual temperature (T), unperturbed, 
    %     ensuring historical correlations with P are preserved*
    % (d) Perturb annual rainfall - change to long term mean
    % (e) Perturb temperature (T) timeseries - change to long term mean
    % (f) Disaggregate P and T timeseries to monthly by method of fragments
    % (g) Convert T to PET using historic regressions
    % (h) Perturb seasonality of rainfall
    % (i) Run rainfall-runoff model to get streamflow
    
    % note, the last two outputs can be considered optional and are as follows: 
    % - TS_SoilMoisture: this output ensures the framework provides
    %   sufficient information to run with the monthly-to-daily disaggregation
    %   technique of John et al. (under review, Rethinking modelling approaches: 
    %   can monthly models replace daily models for hydrological assessments 
    %   under climate change? Submitted to Journal of Hydrology late 2020)
    % - extra: this provides lots of extra information from each step of
    %   the process. 
    
    % (a) Generate timeseries (TS) for high frequency component of annual precipitation, unperturbed
    [TS_P_ann_HighFreq, extra_a] = GetP_Matalas(HistoricData, info);
    
    % (b) Generate low frequency component of annual precipitation, accounting 
    %     for perturbation of Hurst coefficient (low frequency behaviour), and
    %     add the high and low frequency components together
    [TS_P_ann_interim, extra_b] = PerturbP_deltaLowFreqP(TS_P_ann_HighFreq, HistoricData, deltaLowFreqP, info);
    
    % (c) Generate timeseries of annual temperature (T), unperturbed, 
    %     ensuring historical correlations with P are preserved (see notes above) 
    [TS_T_ann_interim, extra_c] = GetTfromP(TS_P_ann_interim, HistoricData, info);
    
    % (d) Perturb long term mean of annual rainfall
    [TS_P_ann, extra_d] = PerturbPmean(TS_P_ann_interim, deltaP, info);    
    
    % (e) Perturb long term mean of temperature
    [TS_T_ann, extra_e] = PerturbT(TS_T_ann_interim, deltaT, info); 
    
    % (f) Disaggregate P and T timeseries to monthly by method of fragments
    [TS_P_interim, TS_T, extra_f] = Disag_AnnToMonthly(TS_P_ann, TS_T_ann, HistoricData, info);
    
    % (g) Convert T to PET using historic regressions
    [TS_PET, extra_g] = GetPETfromT(TS_T, HistoricData, info);
    
    % (h) Perturb seasonality of rainfall
    [TS_P, extra_h] = PerturbP_Seasonality(TS_P_interim, deltaSeasonality, info);
    
    % (i) Get streamflow by running pre-calibrated rainfall runoff model for a set of representative catchments 
    [TS_Q, TS_SoilMoisture, extra_i] = GetQ_RunRRmodel(TS_P, TS_PET, info);
    
    % add ancillary data to data structure 'extra'
    extra = CompileAncillaryData(extra_a, extra_b, extra_c, extra_d, extra_e, extra_f, extra_g, extra_h, extra_i, ... 
                                 TS_P_ann_HighFreq, TS_P_ann_interim, TS_T_ann, TS_T_ann_interim, HistoricData);
    
end

function extra = CompileAncillaryData(extra_a, extra_b, extra_c, extra_d, extra_e, extra_f, extra_g, extra_h, extra_i, TS_P_ann_HighFreq, TS_P_ann_interim, TS_T_ann, TS_T_ann_interim, HistoricData)

    % created 09/10/2019 by Keirnan Fowler, University of Melbourne
    
    extra.extra_a           = extra_a;
    extra.extra_b           = extra_b;
    extra.extra_c           = extra_c;
    extra.extra_d           = extra_d;
    extra.extra_e           = extra_e;
    extra.extra_f           = extra_f;
    extra.extra_g           = extra_g;
    extra.extra_h           = extra_h;
    extra.extra_i           = extra_i;
    extra.TS_P_ann_HighFreq = TS_P_ann_HighFreq;  
    extra.TS_P_ann_interim  = TS_P_ann_interim;   
    extra.TS_T_ann       = TS_T_ann;        
    extra.TS_T_ann_interim  = TS_T_ann_interim;  
    extra.CopyOfHistoricDataInputs = HistoricData;
    
end