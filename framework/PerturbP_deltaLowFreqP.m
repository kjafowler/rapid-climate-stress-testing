
function [TS_P_ann_interim, extra_b] = PerturbP_deltaLowFreqP(TS_P_ann_HighFreq, HistoricData, deltaLowFreqP, info)

    % Created January 2021 by Keirnan Fowler, University of Melbourne   
    
    % This code creates the synthetic annual rainfall timeseries such that 
    % the Hurst coefficient is greater or less than historic, in line with 
    % perturbation indicated by deltaLowFreqP.  A higher Hurst coefficient 
    % means greater persistence in the timeseries which means consecutive 
    % values are more dependant on one another and thus multi-year droughts
    % will tend to be longer, as will multi-year runs of wet conditions.   
    
    % Possible outcomes:
    % deltaLowFreqP = 0 means no change to historic Hurst values; 
    % deltaLowFreqP > 0 means Hurst value goes up by the specified amount;
    % deltaLowFreqP < 0 means Hurst value goes down by the specified amount.
    
    % The 'levers' that this code can 'pull' to acheive the desired Hurst
    % coefficient are the parameters of the model used to generate the low
    % frequency component of the rainfall.  The high frequency component is
    % pre-generated and is passed to this function in TS_P_ann_HighFreq.
    % All Hurst values (whether the context is synethetic data or observed)
    % are calculated on total precipitation rather than the high frequency 
    % or low frequency components separately.   
    
    % The model used to generate the low frequency component is the Broken 
    % Line process.  The two parameters are: 
    % - parameter 1, timescale parameter.
    % - parameter 2, amplitude parameter (standard deviation). 
    % Both parameters must be tuned to give the required response to a given
    % perturbation.  This code assumes a prior routine called 
    % LowFreq_PreAnalysis.m has pre-generated lookup tables of parameters to
    % adopt for a given perturbation.  For the present function, it remains
    % only to apply these parameters to generate the stochastic data.
    
    % for each subarea
    SubareaList = info.SubareaList'; TS_P_ann_LowFreq = []; TS_P_ann_interim = [];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList" 
        
        % get the parameter values to use (with reference to pre-populated tables)
        [param_amplitude, param_timescale] = GetParams(subarea, deltaLowFreqP, info);
        
        % apply the parameter values to generate the low frequency component
        TS_P_ann_LowFreq(:, end+1) = GenTS_BrokenLine(param_timescale, param_amplitude, info.LowFreq_PreAnalysis_Outputs.TS_rand, info.pars);

        % add high and low freqency components together with the long term mean, to get the total rainfall
        TS_P_ann_interim(:, end+1) = P_finalise(TS_P_ann_LowFreq(:, end), TS_P_ann_HighFreq.(subarea{:}), HistoricData.precip_annual.(subarea{:})); 
        
        % remember the value of amplitude and frequency used
        ParamArchive.amplitude.(subarea{:}) = param_amplitude;
        ParamArchive.timescale.(subarea{:}) = param_timescale; 
        
    end
    
    % very rarely, we may get a negative number - check for this
    TS_P_ann_interim(:, :) = max(0.1, TS_P_ann_interim(:, :)); % instead of the negative number, assign the year 0.1 mm of rainfall
    
    % format outputs as tables
    TS_P_ann_LowFreq = array2table(TS_P_ann_LowFreq, 'VariableNames', SubareaList);
    TS_P_ann_interim = array2table(TS_P_ann_interim, 'VariableNames', SubareaList);
    
    extra_b.ParamArchive     = ParamArchive; 
    extra_b.TS_P_ann_LowFreq = TS_P_ann_LowFreq;  
    
end

function [param_amplitude, param_timescale] = GetParams(subarea, deltaLowFreqP, info)
    
    LookupTable = info.LowFreq_PreAnalysis_Outputs.(subarea{:}); 
    
    P = polyfit(LookupTable.perturbation, LookupTable.param_amplitude, 1); param_amplitude = P(1)*deltaLowFreqP + P(2); 
    P = polyfit(LookupTable.perturbation, LookupTable.param_timescale, 1); param_timescale = P(1)*deltaLowFreqP + P(2); 
    
    % note, the reason that we use linear regression not interpolation (or
    % simply reading off the value!) is that random fluctuations when generating
    % the stochastic replicates of the Broken Line process cause irregularities
    % in the relationship, even to the point where it is non-monotonic.  To
    % ensure monotonicity, we take the line of best fit.  
    
end

function CombinedTimeseries = P_finalise(TS_P_ann_LowFreq, TS_P_ann_HighFreq, HistoricData)
    
    LongTermMean = mean(HistoricData); 
    CombinedTimeseries = TS_P_ann_LowFreq ...  % stochastic low frequency component from Broken Line process
                       + TS_P_ann_HighFreq ... % stochastic high frequency component from Matalas
                       + LongTermMean;         % historic long term mean

    % If you are familiar with EMD you will know that EMD also provides a
    % residual which captures the trend in the data (if any).  Often in EMD
    % studies, such trends may be added back to the stochastic data.  
    % However, here it is appropriate that it be removed and not inform the 
    % final stochastic process, because no stress test stochastic timeseries 
    % ever has a trend in it - they are all stationary by design.  Thus, we 
    % are interested only in oscillations, not trends.      
end
