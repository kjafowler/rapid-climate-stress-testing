function [TS_P_ann_HighFreq, extra_a] = GetP_Matalas(HistoricData, info)

    % Created 03/10/2019 by Keirnan Fowler, University of Melbourne
    
    % The purpose of this code is to, for multiple sites, generate a
    % synthetic timeseries of precipitation on an annual timestep using 
    % the Matalas method.  This method preserves the cross-correlation 
    % between sites and the at-site autocorrelation for each individually.  
    
    % reference: 
    % Matalas, N. C. (1967): Mathematical Assessment of Synthetic 
    % Hydrology, Water Resources Research 3(4).  
    
    if ~info.isOctave, rng default; end % reset random seed (not supported in Octave)
    
    UseBoxCox = true; % change to false if you don't want to apply a Box-Cox
                      % transformation to your data - see paper Section ##
    [Precip_BoxCox, BoxCoxParams] = BoxCox(HistoricData, UseBoxCox, info);
    
    % prepare and normalise the data (ie. normalise to a mean of 0 and sd of 1)
    [Matalas_RawData, RawData_stats] = Matalas_prepare(Precip_BoxCox, info); 
    
    % compute Matalas matrices
    [M0, M1, A, B, C] = Matalas_define(Matalas_RawData);
    
    % apply Matalas method to generate stochastic data (in standard normal space)
    Matalas_RawOutput = Matalas_apply(A, B, info.pars.StochRepLen_yrs, info);

    % shift the raw output back to the same space as the input HistoricData.Precip_BoxCox
    Matalas_denormalised = Matalas_denomalise(Matalas_RawOutput, RawData_stats, info); 
    
    % confirm that negative numbers are entirely absent, or you'll get
    % imaginary components when you un-BoxCox.
    testarray = table2array(Matalas_denormalised); testarray_min = min(min(testarray(:, 1:end-1)));
    if testarray_min < 0
        error('Negative numbers present in denormalised Matalas.  Suggest to increase value of "shift" in script kf_BoxCox.m')
    end
    
    % un Box-Cox the data
    TS_P_ann_HighFreq = BoxCox_reverse(Matalas_denormalised, BoxCoxParams, info); 
    
    % add ancillary data to data structure 'extra_a'
    extra_a = CompileAncillaryData(Precip_BoxCox, BoxCoxParams, Matalas_RawData, RawData_stats, M0, M1, A, B, C, Matalas_RawOutput, Matalas_denormalised);

end

function [Matalas_RawData, RawData_stats] = Matalas_prepare(Precip_BoxCox, info)

    % calculate stats for each subarea
    SubareaList = info.SubareaList';  means = []; stds = []; Matalas_RawData_array = [];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in subareas"
        
        % get data and stats for this subarea
        ThisData = Precip_BoxCox.(subarea{:});
        ThisMean = mean(ThisData); 
      % ThisStd  = std(ThisData);    % sample standard deviation (normalised by n-1)
        ThisStd  = std(ThisData, 1); % population standard deviation (normalised by n)
        
        % normalise data
        ThisData_normalised = (ThisData - ThisMean) / ThisStd; 
        
        % append to arrays
        Matalas_RawData_array = [Matalas_RawData_array ThisData_normalised];
        means = [means ThisMean]; 
        stds  = [stds  ThisStd ]; 
    end
    
    % format as tables
    Matalas_RawData    = array2table(Matalas_RawData_array, 'VariableNames', info.SubareaList); 
    RawData_stats.mean = array2table(means,                 'VariableNames', info.SubareaList); 
    RawData_stats.std  = array2table(stds,                  'VariableNames', info.SubareaList); 

end

function Matalas_denormalised = Matalas_denomalise(Matalas_RawOutput, RawData_stats, info)
    
    % get table column numbers
    SubareaList = info.SubareaList'; Matalas_denormalised_array = [];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"
        
        % get data
        NormalisedData = Matalas_RawOutput.(subarea{:});
        ThisMean = RawData_stats.mean.(subarea{:});
        ThisStd  = RawData_stats.std.(subarea{:});
        
        % denormalise
        DenormalisedData = (NormalisedData * ThisStd) + ThisMean; 
        
        % append to array
        Matalas_denormalised_array = [Matalas_denormalised_array DenormalisedData];
        
    end
    
    % convert array to table
    Matalas_denormalised = array2table(Matalas_denormalised_array, 'VariableNames', info.SubareaList);

end

function extra_a = CompileAncillaryData(Precip_BoxCox, BoxCoxParams, Matalas_RawData, RawData_stats, M0, M1, A, B, C, Matalas_RawOutput, Matalas_denormalised)
    extra_a.Precip_BoxCox = Precip_BoxCox;
    extra_a.BoxCoxParams = BoxCoxParams;
    extra_a.Matalas_RawData = Matalas_RawData;
    extra_a.RawData_stats = RawData_stats;
    extra_a.M0 = M0;
    extra_a.M1 = M1;
    extra_a.A = A;
    extra_a.B = B;
    extra_a.C = C;
    extra_a.Matalas_RawOutput = Matalas_RawOutput;
    extra_a.Matalas_denormalised = Matalas_denormalised;
end

