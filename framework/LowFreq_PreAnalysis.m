function [LowFreq_PreAnalysis_Outputs, extra] = LowFreq_PreAnalysis(HistoricData, info)

    % created February 2021 by Keirnan Fowler, University of Melbourne
    
    % NOTES ABOUT STOCHASTIC GENERATION FOR THE BASE CASE (NO PERTURBATION)
    % 1.  This framework explicity separates the annual precipitation into
    %     high frequency and low frequency components, and stochastically 
    %     generates each component separately.  
    % 2.  The model used for the low frequency component is called the
    %     'Broken Line Process'.  It has two parameters that need to be 
    %     decided: the timescale parameter and the amplitude parameter.  
    % 3.  To ensure all subareas synchronously go into and out of multiyear 
    %     wet and dry periods, all subareas need to have the same timescale
    %     parameter.  
    % 4.  For the base case, the timescale parameter is decided based on
    %     the number of zero crossings in the combined historic low
    %     frequency series (ie. combined across subareas).  The formula
    %     used is in Supplementary Material ##.  
    % 5.  The amplitude parameter is then decided for each subarea
    %     individually.  The basis for deciding is that the Historic Hurst
    %     exponent value for that subarea must be preserved in the final 
    %     stochastic data (high+low combined).  
    
    
    % NOTES ABOUT PERTURBATION
    % 6.  As part of the stress test, the low frequency component can be
    %     perturbed.  
    % 7.  Perturbations in low frequency component are quantified as a
    %     change in Hurst exponent.  
    % 8.  The change in Hurst exponent is acheived solely by changing the
    %     parameters of the Broken Line process - no change to high
    %     frequency stochastic generation is made. 
    % 9.  An increase in Hurst can be achieved by changing the low frequency
    %     component so that it has a longer period, or higher amplitude
    %     oscillations, or both.  Here, we opt for both.  
    % 10. To understand how this is acheived, imagine the following (for a
    %     given subarea): a 2D plot, x = timescale parameter; y = amplitude 
    %     parameter.  The base case parameter
    %     values from 4. and 5. plot as a point in this space.  For a given
    %     Hurst value, we can connect all the points in the space that
    %     result in that value (post being recombined with the high 
    %     frequency component).  For a higher Hurst value, the isoline will
    %     be up and to the right (longer period and/or higher amplitude).
    % 11. During perturbation, the parameters are changed so that the
    %     direction of change is perpendicular to the isolines (or as
    %     close to this as possible given the restrictions in 12).  As per 
    %     10, this means the direction will be up and to the right, or put
    %     another way it means the angle will be between 0 and 90 degrees.
    % 12. As per 3, all regions need to have the same timescale parameter.
    %     The value of this parameter will change with perturbation of 
    %     Hurst, but the key point is that for a given perturbation all 
    %     subareas need to have the same timescale parameter value. As per
    %     5, this does not apply to the amplitude parameter.  In fact, the
    %     amplitude parameter can change to ensure the desired (perturbed)
    %     Hurst value is obtained in each subarea.  
    
    % THIS CODE
    % This code is intended as a pre-analysis to be undertaken prior to
    % stochastic data generation.  It is faster to do the routines herein
    % once only at the start, rather than repeating them for every possible
    % perturbation.  
    
    % OUTPUT OF THIS CODE
    % Main output is a set of tables (one for each subarea) where the columns are:
    % - perturbation
    % - resulting Hurst value
    % - timescale parameter value
    % - amplitude parameter value
    % The idea is to trial many different perturbation values and then at
    % run time the required parameter values can be read off the table, or 
    % interpolated if the exact perturbation value has not been tested.   
    
    % STEPS UNDERTAKEN IN THIS CODE
    % (a) For each subarea, calculate the historic Hurst value.
    % (b) Calculate the base case timescale parameter as per 4.  
    % (c) For each subarea, stochastically generate the high frequency 
    %     component in exactly the same manner as happens at run time.
    % (d) For each subarea, calculate the base case amplitude parameter as per 5.  
    % (e) For each subarea, test every angle between 1 degree and 89 degrees 
    %     (see point 11), in increments of 1 degree.  The test involves 
    %     moving a set distance (pars.test_distance) in the specified 
    %     direction and recording the change in Hurst.  
    % (f) Decide on a direction (to be applied across all subareas) that 
    %     gives the greatest change in Hurst value, on average, across all 
    %     the subareas.  
    % (g) Calculate the gradient [?Hurst / dist. in 2D space] for each
    %     subarea separately and then take the average across all subareas.
    % (h) Use the gradient to populate the third column of the tables (ie.
    %     for each perturbation, assign a timescale parameter)
    % (i) For each perturbation, given the timescale parameter, solve for
    %     the exact value of the amplitude parameter that will give the
    %     desired (perturbed) Hurst value.  Do this separately for each
    %     subarea.  
    
    % set perturbations to test (defined as change in Hurst value)
    % These don't neccesarily have to correspond to the adopted perturbations - 
    % if not, interpolation is used.
    output_table.perturbation = [-0.03 -0.015 0 +0.015 +0.03 +0.045 +0.06 +0.075];
    
    %% (a) For each subarea, calculate the historic Hurst value.
    HistoricHurst = GetHistoricHurst(HistoricData, info);
    
    %% (b) Calculate the base case timescale parameter as per 4.  
    param_timescale_hist = GetParamTimescaleHist(HistoricData, info);  
    
    %% (c) For each subarea, stochastically generate the high frequency 
    %     component in exactly the same manner as happens at run time.
    TS_P_ann_HighFreq = GetP_Matalas(HistoricData, info);
    
    %% (d) For each subarea, calculate the base case amplitude parameter as per 5.  
    disp('Starting base case parameter calculations...')
    [param_amplitude_hist, TS_rand] = GetParamAmp_unpert(HistoricData, param_timescale_hist, info, TS_P_ann_HighFreq, HistoricHurst);
    
    %% (e) For each subarea, test every angle between 1 degree and 89 degrees 
    %     (see point 11), in increments of 1 degree.  The test involves 
    %     moving a set distance (pars.test_distance) in the specified 
    %     direction and recording the change in Hurst.  
    disp('Done. Starting directional tests for Broken Line parameters...')
    info.pars.test_distance = 0.1; % one tenth of a unit
    DirectionalTestResults = DirectionalTest(HistoricData, param_timescale_hist, info, TS_P_ann_HighFreq, HistoricHurst, param_amplitude_hist, TS_rand); 
    
    %% (f) Decide on a direction (to be applied across all subareas) that 
    %     gives the greatest change in Hurst value, on average, across all 
    %     the subareas.  
    % and
    %  (g) Calculate the gradient [?Hurst / dist. in 2D space] for each
    %     subarea separately and then take the average across all subareas.
    [upslope_angle_deg, grad] = FinaliseDirection(DirectionalTestResults, info); 
    
    %% (h) Use the gradient and angle to populate the second column of the tables (ie.
    %     for each perturbation, assign a timescale parameter)
    output_table.timescale = GetTimescale(output_table.perturbation, grad, upslope_angle_deg, param_timescale_hist); 
    
    %% (i) For each perturbation, given the timescale parameter, solve for
    %     the exact value of the amplitude parameter that will give the
    %     desired (perturbed) Hurst value.  Do this separately for each
    %     subarea.  
    disp('Done.  Finalising Broken Line parameters...')
    [output_table.amplitude, output_table.Hurst] = GetAmplitude(HistoricData, info, TS_P_ann_HighFreq, HistoricHurst, TS_rand, output_table.timescale, output_table.perturbation); 
    
    %% create final output
    LowFreq_PreAnalysis_Outputs = CreateFinalOutput(output_table, TS_rand, info);
    extra = CreateExtraOutput(info, DirectionalTestResults, upslope_angle_deg, grad); 

end

%% functions

function HistoricHurst = GetHistoricHurst(HistoricData, info)
    SubareaList = info.SubareaList'; HistoricHurst = [];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"
        ThisHistHurst = GetHurst(HistoricData.precip_annual.(subarea{:}), false);  
        HistoricHurst(end+1) = ThisHistHurst;
    end
    
    % save as table
    subarea = SubareaList'; HistoricHurst = HistoricHurst'; 
    HistoricHurst = table(subarea, HistoricHurst); 
end

function param_timescale_hist = GetParamTimescaleHist(HistoricData, info) 
    TimeseriesLength = size(HistoricData.precip_ann_LowFreq.adopted_combined, 1); 
    NumZeroCrossings = GetNumZeroCrossings(HistoricData.precip_ann_LowFreq.adopted_combined);
    param_timescale_hist = TimeseriesLength / (2*NumZeroCrossings + 1); % see Supplementary Material ##
end

function [param_amplitude_hist, TS_rand] = GetParamAmp_unpert(HistoricData, param_timescale_hist, info, TS_P_ann_HighFreq, HistoricHurst);
    
    % generate TS_rand, a single timeseries of random deviations for use across all regions. 
    n = size(TS_P_ann_HighFreq, 1); TS_rand = rand(n+1, 1); 
    
    SubareaList = info.SubareaList'; param_amplitude_hist = [];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"    
        
        pos = strcmp(HistoricHurst.subarea, subarea{:}); 
        ThisHistHurst = HistoricHurst.HistoricHurst(pos);
        
        % do a search to find the best param_amplitude_hist so that the
        % historic Hurst value is matched    
        fun = @(param_amplitude) abs(ThisHistHurst - GetStochHurstGivenPars(param_timescale_hist, param_amplitude, TS_rand, info, TS_P_ann_HighFreq, HistoricData, subarea)); 
        initial_val = 100; 
        options = optimset('TolFun',1e-3, 'TolX', 1e-2);
        param_amplitude_hist(end+1) = fminsearch(fun, initial_val, options); 
    
    end
    
    % save as table
    subarea = SubareaList'; param_amplitude_hist = param_amplitude_hist'; 
    param_amplitude_hist = table(subarea, param_amplitude_hist);     
    
end

function DirectionalTestResults = DirectionalTest(HistoricData, param_timescale_hist, info, TS_P_ann_HighFreq, HistoricHurst, param_amplitude_hist, TS_rand)
    
    % Notes: 
    % The objective is to find the angle (direction in parameter space) that
    % gives the greatest increase in Hurst.  While fminsearch would have
    % been a good option here, the objective function surface is often
    % irregular** and thus may be multi-modal, so we use a more manual 
    % approach:  We manually test every angle between 1 and 89 degrees, in 
    % 1 degree increments.  
    
    % ** I think because higher values of param_timescale_hist incorporate 
    %    progressively less Broken Line segments in the stochastic timeseries, 
    %    and the removal of each segment has a unique affect on Hurst.      
    
    SubareaList = info.SubareaList'; 
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"    
        
        pos = strcmp(HistoricHurst.subarea, subarea{:}); 
        ThisHistHurst = HistoricHurst.HistoricHurst(pos);
        this_param_amplitude_hist = param_amplitude_hist.param_amplitude_hist(pos); 
        
        AnglesToTest = 1:89; % anywhere between 1 degree and 89 degrees, in one-degree increments
        for i = 1:size(AnglesToTest, 2)
            angle = AnglesToTest(i); 
            [HurstChange(i), NewHurst(i)] = GetHurstChangeForSlope(angle, param_timescale_hist, this_param_amplitude_hist, ThisHistHurst, TS_rand, info, TS_P_ann_HighFreq, HistoricData, subarea); 
        end
        DirectionalTestResults.HurstChange.(subarea{:}) = HurstChange; 
        DirectionalTestResults.NewHurst.(subarea{:}) = NewHurst; 
        disp(['Done for subarea ' subarea{:}])
    end

end

function [upslope_angle_deg, grad] = FinaliseDirection(DirectionalTestResults, info)

    plotting = false; % make true if you want to see a plot
 
    SubareaList = info.SubareaList'; array_AllSubareas = [];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"
        array_AllSubareas(1:89, end+1) = DirectionalTestResults.HurstChange.(subarea{:});
        if plotting, plot(DirectionalTestResults.HurstChange.(subarea{:}), '-o'); hold on; end
    end
    mean_AllSubareas = mean(array_AllSubareas, 2);
    [delta_Hurst, upslope_angle_deg] = max(mean_AllSubareas);
    if plotting, plot(mean_AllSubareas, 'LineWidth', 5); end
    
    % delta_Hurst is the average increase in Hurst over pars.test_distance.  
    % However, we need the gradient, which is the change over a single unit of distance.  
    grad = delta_Hurst / info.pars.test_distance;
    
end

function timescale_param = GetTimescale(list_of_perturbations, grad, upslope_angle_deg, param_timescale_hist)
    
    % Use the gradient to assign a timescale parameter for each perturbation.  
    
    % Logic is as follows:
    % if grad is 0.2 it means the Hurst changes by 0.2 for each unit
    % of distance traversed in 2D parameter space.  So if we need a
    % perturbation of 0.3 we need to traverse 1.5 units in 2D space.
    % However, the timescale parameter indicates distance in the x direction, 
    % whereas the above distance is travelled in the direction indicated by
    % upslope_angle_deg.  Thus we need to multiply by cos(angle) to get the
    % distance on the x-axis.  
    
    NumPert = size(list_of_perturbations, 2); 
    for iPert = 1:NumPert
        ThisPert = list_of_perturbations(iPert);
        DistTravelled_angle = ThisPert / grad;
        DistTravelled_xdir = DistTravelled_angle * cos(upslope_angle_deg / (180/pi())); 
        timescale_param(iPert) = param_timescale_hist + DistTravelled_xdir;
    end
    
end

function [amplitude, Hurst] = GetAmplitude(HistoricData, info, TS_P_ann_HighFreq, HistoricHurst, TS_rand, timescales, perturbations)
    
    SubareaList = info.SubareaList';  
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"    

        % initialisation
        pos = strcmp(HistoricHurst.subarea, subarea{:}); 
        ThisHistHurst = HistoricHurst.HistoricHurst(pos);
        amplitude.(subarea{:}) = []; Hurst.(subarea{:}) = []; 
        
        NumPert = size(perturbations, 2);
        for iPert = 1:NumPert
            
            % key values for input into search algorithm
            ThisPert = perturbations(iPert); 
            TargetHurst = ThisHistHurst + ThisPert; 
            ThisParamTimescale = timescales(iPert); 
            
            % given the above, we now conduct a search to identify the
            % value of the amplitude parameter that gives us TargetHurst. 
            fun = @(param_amplitude) abs(TargetHurst - GetStochHurstGivenPars(ThisParamTimescale, param_amplitude, TS_rand, info, TS_P_ann_HighFreq, HistoricData, subarea)); 
            initial_val = 50; 
            options = optimset('TolFun',1e-3, 'TolX', 1e-2);
            ThisParamAmplitude = fminsearch(fun, initial_val, options); 
            
            % store information in function output tables
            amplitude.(subarea{:})(end+1) = ThisParamAmplitude; 
            Hurst.    (subarea{:})(end+1) = TargetHurst; 
            
        end
        disp(['Done for subarea ' subarea{:}])
    end
end

function LowFreq_PreAnalysis_Outputs = CreateFinalOutput(output_table, TS_rand, info)

    SubareaList = info.SubareaList'; 
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"
        perturbation    = output_table.perturbation';           
        Hurst_NewVal    = output_table.Hurst.(subarea{:})';     
        param_timescale = output_table.timescale';              
        param_amplitude = output_table.amplitude.(subarea{:})'; 
        
        LowFreq_PreAnalysis_Outputs.(subarea{:}) = table(perturbation, Hurst_NewVal, param_timescale, param_amplitude); 
    end
    LowFreq_PreAnalysis_Outputs.TS_rand = TS_rand;
    
end

function extra = CreateExtraOutput(info, DirectionalTestResults, upslope_angle_deg, grad)
    extra.test_distance = info.pars.test_distance;
    extra.DirectionalTestResults = DirectionalTestResults;
    extra.upslope_angle_deg = upslope_angle_deg;
    extra.grad = grad;  
end

function [HurstChange, NewHurst] = GetHurstChangeForSlope(angle, param_timescale_hist, param_amplitude_hist, HistHurst, TS_rand, info, TS_P_ann_HighFreq, HistoricData, subarea)
    
    % determine the change in Hurst if we traverse a set distance** in the
    % direction specified.  This code is best read by visualing the 2D
    % space with the timescale parameter on the x-axis and the amplitude
    % parameter * 0.01 on the y-axis.  
    
    % ** for the set distance we use one tenth of a unit
    
    param_timescale_test =       param_timescale_hist      + info.pars.test_distance*cos(angle / (180/pi())); 
    param_amplitude_test = 100*((param_amplitude_hist/100) + info.pars.test_distance*sin(angle / (180/pi()))); 
    
    NewHurst = GetStochHurstGivenPars(param_timescale_test, param_amplitude_test, TS_rand, info, TS_P_ann_HighFreq, HistoricData, subarea);
    
    HurstChange = NewHurst - HistHurst; 
    
end

function NumCrossings = GetNumZeroCrossings(ts)
    NumCrossings=0; 
    for k=2:size(ts,1) 
        if ts(k) * ts(k-1)<0
            NumCrossings=NumCrossings+1;
        end
    end
end

function [StochHurst, TS_BrokenLine] = GetStochHurstGivenPars(param_timescale, param_amplitude, TS_rand, info, TS_P_ann_HighFreq, HistoricData, subarea)
    
    % generate the stochastic timeseries for the low frequency component, 
    % using the Broken Line process with the parameter values specified
    TS_BrokenLine = GenTS_BrokenLine(param_timescale, param_amplitude, TS_rand, info.pars);
    
    % sum the stochastic low and high frequency components together, along with the mean annual precipitation 
    TotalP = TS_BrokenLine' + TS_P_ann_HighFreq.(subarea{:}) + mean(HistoricData.precip_annual.(subarea{:})); 
    
    % calculate the Hurst value of the combined timeseries
    % StochHurst = GetHurstMean(TotalP, size(HistoricData.precip_annual, 1)); 
    StochHurst = GetHurstMean(TotalP, size(HistoricData.precip_annual, 1)); 
end

function HurstCoeff = GetHurstMean(data, HistSeriesLength)

    % get an average Hurst Coefficient for an arbitrarily long period, by
    % splitting it into shorter segments of equal length to the historic
    % data and calculating the Hurst Coefficient separately for each of the
    % segments, then taking the average.  

    
    len = size(data, 1);
    
    ind(1) = 1; this_count = 1;
    for i = 2:len
        if this_count == HistSeriesLength % divide into segments of length HistSeriesLength
            ind(i) = ind(i-1)+1; this_count = 1;
        else
            ind(i) = ind(i-1); this_count = this_count + 1; 
        end
    end
    
    % calculate the number of segments
    NumSegments = max(ind); 
    
    % calculate the Hurst coefficient for each of the segments separately
    for iSegment = 1:NumSegments
        ind2 = ind==iSegment;
        Hurst_all(iSegment) = GetHurst(data(ind2), 'PlottingOff');
    end
    
    HurstCoeff = mean(Hurst_all); 

end

function HurstCoeff = GetHurst(data, PlotToggle)

    % created 16/02/2020 by Keirnan Fowler, University of Melbourne
    
    % calculate the original version of the Hurst Coefficient (ie. the one
    % hydrologists use, not economists).  Only input is a timeseries of
    % data.  
    
    % divide into partial periods
    [NumDurations, durations, NumPartials, ind] = GetPartials(data);
    
    % calculate rescaled range for each partial period
    RSranges = GetRescaledRanges(data, ind, NumDurations, NumPartials); 
    
    % calculate Hurst exponent by fitting power law to results
    % and, optionally, produce plots
    HurstCoeff = GetHurstAndPlot(RSranges, durations, NumDurations, NumPartials, PlotToggle);
    
end

function [NumDurations, durations, NumPartials, ind] = GetPartials(data)
    
    n = size(data, 1);
    
    % Assemble information about segment duration when the full period is
    % divided into N, N/2, N/4, N/8 etc...

    NextDuration = n; NumDivs = 1; iDuration = 0;
    while NextDuration>=2
        
        % information about this duration
        iDuration = iDuration + 1; 
        durations(iDuration) = NextDuration;
        NumPartials(iDuration) = floor(n/durations(iDuration));
        
        % information about the next shortest duration (we will not proceed
        % unless there are sufficient years)
        NumDivs = NumDivs*2; 
        NextDuration = floor(n/NumDivs);
    end
    NumDurations = iDuration;
    
    
    % create an array (ind) that indicates the years that belong to each segment
    ind = zeros(NumDurations, n);
    for iDuration = 1:NumDurations
        pos_end= 0;
        for iPartial = 1:NumPartials(iDuration)
            pos_start = pos_end + 1; 
            pos_end = pos_start + durations(iDuration) - 1; 
            ind(iDuration, pos_start:pos_end) = iPartial; 
        end        
    end    
end

function RSranges = GetRescaledRanges(data, ind, NumDurations, NumPartials)
    n = size(data, 1); 
    for iDuration = 1:NumDurations
        for iPartial = 1:NumPartials(iDuration)
            
            % get this partial timeseries
            ThisInd = find(ind(iDuration, :)==iPartial);
            ThisTS = data(ThisInd); 
            
            % calculate the mean and standard deviation
            ThisMean = mean(ThisTS); ThisSD = std(ThisTS);
            
            % create a mean-adjusted series
            ThisTS = ThisTS - ThisMean; 
            
            % calculate the cumulative deviate series
            ThisTS = cumsum(ThisTS);
            
            % compute the range
            ThisRange = range(ThisTS);

            % compute the rescaled range
            RSranges(iDuration, iPartial) = ThisRange / ThisSD;
        end
    end
end

function HurstCoeff = GetHurstAndPlot(RSranges, durations, NumDurations, NumPartials, PlotToggle)

    %% Get Hurst exponent
    
    % get relevant ordinates in log-log space
    for iDuration = 1:NumDurations
        x(iDuration) = log(durations(iDuration));
        y(iDuration) = log(mean(RSranges(iDuration, 1:NumPartials(iDuration))));
    end
    
    % line of best fit
    [m, c] = LinReg_get_m_and_c(x, y);
    
    % Hurst exponent is slope of log-log line
    HurstCoeff = m;

    %% optional - do plotting
    DoPlots = false;
    if strcmp(PlotToggle, 'PlottingOn'), DoPlots = true; end
    
    if DoPlots
        
        % remember the x and y values we just calculated
        x_log = x; y_log = y;
        
        % initiate figure
        figure; 
        
        % subplot 1: rescaled ranges versus partial record length
        subplot(1, 2, 1); hold on;

        % individual values for each instance of a partial record
        for iDuration = 1:NumDurations
            x = repmat(durations(iDuration), NumPartials(iDuration), 1);
            y = RSranges(iDuration, 1:NumPartials(iDuration));
            scatter(x, y); 
        end        

        % average values for each value of record length (unlogged)
        x = []; y = [];
        for iDuration = 1:NumDurations
            x(iDuration) = durations(iDuration);
            y(iDuration) = mean(RSranges(iDuration, 1:NumPartials(iDuration)));
        end
        scatter(x, y, 30, 'k', 'filled'); 
        grid on; xlabel('length of partial record (years)'); ylabel('rescaled range (-)');

        % subplot 2: averages only, in log-log space, with fitted line
        subplot(1, 2, 2); hold on;
        scatter(x_log, y_log, 30, 'k', 'filled'); 
        grid on; xlabel('log_{n}(length in years)'); ylabel('log_{n}(rescaled range)')

        % line of best fit
        plot([1.1 5], [m*1.1+c, m*5+c], '-r');
        
        % labelling
        text(3, 1, ['Hurst coefficient = ' num2str(HurstCoeff)]);
    end        
    
end

function [m, c] = LinReg_get_m_and_c(x, y)
    P = polyfit(x,y,1);
    m = P(1); c = P(2);
end
