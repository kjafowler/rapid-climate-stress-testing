
function example_plotting(StochPertData, HistoricData, info)

    % Created January 2021 by Keirnan Fowler, University of Melbourne, fowler.k@unimelb.edu.au

    % Integrated framework for rapid climate stress testing on a monthly timestep
    % by Keirnan Fowler, Natasha Ballis, Avril Horne, Andrew John, Rory Nathan and Murray Peel 

    % Licence: CC BY 3.0 - see https://creativecommons.org/licenses/by/3.0/au/

    % The purpose of this function is to produce the plots for the example
    % given in the paper.  Note: in the more complex figures (eg. Figure 6a,  
    % 10, 12) this code provides only the base data that is plotted 
    % and *not* extras such as figure labels, axis labels, and legends.  
    % This is to keep the code simple and also because the authors added 
    % these features in postprocessing using Inkscape (https://inkscape.org/).  
    
    if ~info.isOctave

        % Figure 7a: separation of historic annual rainfall into low and high frequency components 
        Figure7a(HistoricData, info);
        disp('Figure 7a done.'); 

        % Figure 7b: combination of low frequency component across subareas
        Figure7b(HistoricData, info); 
        disp('Figure 7b done.'); 

        % Figure 8: seasonality of precipitation, temperature, PET and streamflow
        Figure8(HistoricData, StochPertData);
        disp('Figure 8 done.'); 

        % Figure 9: Match of stochastic baseline streamflow with historic for
        % annual and 5 year totals
        Figure9(HistoricData, StochPertData, info)
        disp('Figure 9 done.'); 

        % Figure 11: Evaluation of stochastic data
        HistRecordLength = size(HistoricData.flow_annual.ReachInflows_ML, 1);
        Figure11(StochPertData, info, HistRecordLength);
        disp('Figure 11 done.');

        % Figure 12: Stress testing results
        Figure12(); % note: no inputs, because this figure uses dummy data
        disp('Figure 12 done')

        % optionally, produce the quality evaluation plots in the supplementary material
        ProduceSuppMattPlots = true;
        if ProduceSuppMattPlots, SuppMattPlots(HistoricData, StochPertData, info); disp('Supplementary plots done.'); end

    else
        warning('To users of Octave: plots for the paper were produced using Matlab code which has not been tested in Octave.  Please write your own plotting routines as required.')
    end

end

function Figure7a(HistoricData, info)
    
    figure

    % for each subarea
    SubareaList = info.SubareaList'; iSubarea = 0; 
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"    
             
        iSubarea = iSubarea + 1; 
        
        % plotting - original signal
        x = HistoricData.precip_annual.Year';
        OriginalSignal = HistoricData.precip_annual.(subarea{:});
        subplot(info.pars.NumSubareas, 3, 1+(iSubarea-1)*3); hold on;
        plot(x, OriginalSignal, 'k'); ylim([0 2000])
        set(gca, 'XTickLabels', [], 'YTickLabels', []);  
        
        % plotting - remove Intrinsic Mode Functions, in two steps
        BeforeRemoval = OriginalSignal; 
        for iPlot = 2:3
            
            switch iPlot
                case 2, component_to_remove = HistoricData.precip_ann_HighFreq.(subarea{:});
                case 3, component_to_remove = HistoricData.precip_ann_LowFreq. (subarea{:});
            end

            % subtract the component from the starting series
            AfterRemoval = BeforeRemoval - component_to_remove; 

            % prepare matrices to plot area in between curves
            x2 = [x, fliplr(x)]; inBetween = [AfterRemoval; flipud(BeforeRemoval)]; 

            % plotting
            subplot(7, 3, iPlot+(iSubarea-1)*3); hold on;
            fill(x2, inBetween, [1.0 0.4 0.4], 'LineStyle','none');
            plot(x, AfterRemoval, 'k'); 
            ylim([0 2000]); % user to specify desired y limits
            set(gca, 'XTickLabels', [], 'YTickLabels', []);  

            % prepare for next IMF
            BeforeRemoval = AfterRemoval; 
        end
    end
    
    % print the figure
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 15 17]); 
    print('-painters','-depsc', 'out\figures\Figure7a_raw'); saveas(gcf, 'out\figures\Figure7a_raw', 'png');   
    close
    
end

function Figure7b(HistoricData, info)

    % initialise some vars
    SubareaList = info.SubareaList'; 
    LineStyles = {':', '--', '-', '-.', '--', '-.', ':'};
    x = HistoricData.precip_annual.Year; 
    iSubarea = 0; 
    figure; hold on; 
    
    % line for each subarea
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"    
        y = HistoricData.precip_ann_LowFreq.(subarea{:});
        iSubarea = iSubarea + 1; 
        plot(x, y, LineStyles{iSubarea}, 'LineWidth', 1.5); 
    end
    
    % line for weighted average across subareas
    y = HistoricData.precip_ann_LowFreq.adopted_combined; 
    plot(x, y, 'k', 'LineWidth', 3); 
    
    % plot formatting
    ylabel('Low frequency component (mm)');
    legend('sub-area A', 'sub-area B', 'sub-area C', 'sub-area D', 'sub-area E', 'sub-area F', 'sub-area G', 'Adopted', 'location', 'southoutside'); 
    xlim([1910 2020]); set(gca, 'xTick', 1920:20:2020); 
    ylim([-300 300]); set(gca, 'yTick', -300:100:300); grid on
    
    % save the figure
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 10 15]); 
    print('-painters','-depsc', 'out\figures\Figure7b'); saveas(gcf,'out\figures\Figure7b', 'png'); 
    close

end

function Figure8(HistoricData, StochPertData)
    
    iSubarea = 2;  % subarea B (but this can be changed so as to focus on a 
                   % different subarea as required)

    figure; 

    for iDataType = 1:4

        switch iDataType
            case 1, data.hist = HistoricData.precip_monthly;                           data.stoch = StochPertData.BaseCase.TS_P;                plottitle = '(a) precip.';    yaxislabel = 'precipitation (% in each month)'; 
            case 2, data.hist = HistoricData.T_monthly;                                data.stoch = StochPertData.BaseCase.TS_T;                plottitle = '(b) Tmax.';      yaxislabel = 'Tmax (°C, average monthly)'; 
            case 3, data.hist = HistoricData.PET_monthly;                              data.stoch = StochPertData.BaseCase.TS_PET;              plottitle = '(c) PET';        yaxislabel = 'PET (% in each month)'; 
            case 4, data.hist = HistoricData.flow_monthly.RepresentativeCatchments_mm; data.stoch = StochPertData.BaseCase.TS_Q_ReprCatchments; plottitle = '(d) streamflow'; yaxislabel = 'streamflow (% in each month)';
        end

        % for all except T, convert to monthly %
        if iDataType ~= 2
            data.hist  = ConvertToMonthlyPercent(data.hist);
            data.stoch = ConvertToMonthlyPercent(data.stoch);
        end            

        NumCols = 0; bp_data = [];
        for iMonth = 1:12

            % develop index identifying data for this month
            ind_hist  = data.hist.Month == iMonth;  NumMon_hist  = sum(ind_hist); 
            ind_stoch = data.stoch.Month == iMonth; NumMon_stoch = sum(ind_stoch); 

            % get data for this month
            hist_ThisMon  = data.hist {ind_hist,  iSubarea+2}; 
            stoch_ThisMon = data.stoch{ind_stoch, iSubarea+2}; 


            % create nan-filled columns for data to go into
            bp_data (1:NumMon_stoch, (NumCols+1):(NumCols+3)) = nan; 

            % add data to boxplot table
            bp_data (1:NumMon_hist,  NumCols+1) = hist_ThisMon; 
            bp_data (1:NumMon_stoch, NumCols+2) = stoch_ThisMon; 
            bp_ColourGroup    ((NumCols+1):(NumCols+3)) = [2 3 1];

            NumCols = NumCols + 3; 

        end

        colours_array = [  0   0   0; ... % blank series - colour irrelevant
                         0.3 0.8 0.3; ... % historic - green
                         0.5 0.5 0.5];    % stochastic baseline - grey


        subplot(1, 4, iDataType); hold on;
        bp_labels = {'J', '', '', 'F', '', '','M', '', '','A', '', '','M', '', '','J', '', '','J', '', '','A', '', '','S', '', '','O', '', '','N', '', '','D', '', ''}; 
        bp = boxplot(bp_data, 'PlotStyle', 'compact', 'ColorGroup', bp_ColourGroup, 'Labels', bp_labels, 'Colors', colours_array, 'Symbol', '');

        % bp = boxplot(bp_data, 'PlotStyle', 'compact', 'ColorGroup', bp_ColourGroup, 'Labels', bp_labels, 'Colors', colours_array, 'Symbol', '');

        % fiddle with whisker length - make them 5th and 95th
        WhiskPct_lwr = 5; WhiskPct_upp = 95; 
        lines = get(gca, 'children');
        uw = findobj(lines, 'tag', 'Whisker'); numBx = size(uw, 1); assert(numBx == size(bp_data, 2));
        for i = 1:numBx, uw(i).YData = [prctile(bp_data(:, numBx-i+1), WhiskPct_lwr) prctile(bp_data(:, numBx-i+1), WhiskPct_upp)]; end            

        % other formatting
        ylabel(yaxislabel)
        title(plottitle); 
        switch iDataType
            case 1, ylim([0 20]);
            case 2, ylim([0 30]);
            case 3, ylim([0 20]);
            case 4, ylim([0 30]);
        end

        if iDataType == 4
            % to trick Matlab into creating a suitable legend, draw
            % two additional lines out of site
            yl = ylim; y_OutOfSight = yl(2) + 50; 
            plot([10 20], [y_OutOfSight y_OutOfSight], 'Color', [0.3 0.8 0.3], 'LineWidth', 5); 
            plot([10 20], [y_OutOfSight y_OutOfSight], 'Color', [0.5 0.5 0.5], 'LineWidth', 5); 

            % now add legend
            legend('historic', 'stochastic');
        end        

    end

    % save the figure
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 28 15]); 
    print('-painters','-depsc', 'out\figures\Figure8'); saveas(gcf,'out\figures\Figure8', 'png'); 
    close
    
end

function Figure9(HistoricData, StochPertData, info)
    
    % aggregate stochastic data up to annual
    StochPertData.BaseCase = AggregateToAnnual_stoch(StochPertData.BaseCase, info.pars);
    
    % specify data sources
    hist = HistoricData.flow_annual.ReachInflows_ML.Reach1 ...
         + HistoricData.flow_annual.ReachInflows_ML.Reach2 ...
         + HistoricData.flow_annual.ReachInflows_ML.Reach3 ...
         + HistoricData.flow_annual.ReachInflows_ML.Reach4; % sum of all unreg reaches
    
    stoch = StochPertData.BaseCase.TS_annual_Q_ReachInflows.Reach1 ...
          + StochPertData.BaseCase.TS_annual_Q_ReachInflows.Reach2 ...
          + StochPertData.BaseCase.TS_annual_Q_ReachInflows.Reach3 ...
          + StochPertData.BaseCase.TS_annual_Q_ReachInflows.Reach4; 
    
    % format stoch data into replicates of same length as historic record
    HistRecordLength = size(HistoricData.flow_annual.ReachInflows_ML, 1);
    stoch_rearranged = RearrangeData(HistRecordLength, stoch); 
    NumReplicates = size(stoch_rearranged, 2); 
    
    % settings
    overlapping = true; 
    ReqdDurations = [1 5]; % require separate plots for 1 and 5 year runs
    NumDur = size(ReqdDurations, 2); 
    
    
    for iDur = 1:NumDur
    
        %% Data extraction
        
        % initialise
        ThisDur = ReqdDurations(iDur);
        
        % collate frequency curves
        sorted = true;
        HistRuns = CollateFreqCurve(hist, ThisDur, overlapping, sorted); 
        for iRep = 1:NumReplicates
            StochRuns(:, iRep) = CollateFreqCurve(stoch_rearranged(:, iRep), ThisDur, overlapping, sorted); 
        end

        % extract various percentiles out of the stochastic runs
        percentiles = [1 5 10 20 50 70 90 95 99]; StochPctiles = [];
        for iPct = 1:size(percentiles, 2)
            NumOrd = size(StochRuns, 1); 
            for iOrd = 1:NumOrd
                StochPctiles(iOrd, iPct) = prctile(StochRuns(iOrd, :), percentiles(iPct)); 
            end
        end
        
        %% plotting
        
        subplot(1, 2, iDur); hold on;
        
        % plot historic data
        NumOrds = size(HistRuns, 2); 
        x = 100 * ((1:NumOrds) / NumOrds); x = fliplr(x); 
        scatter(x, HistRuns, 22, 'r', 'filled');    
        
        % plot stochastic percentiles
        plot(x, StochPctiles(:, 2), '--', 'Color', [0.3 0.3 0.3]); % 5th
        plot(x, StochPctiles(:, 5), '-', 'Color', [0.3 0.3 0.3]);  % median
        plot(x, StochPctiles(:, 8), '--', 'Color', [0.3 0.3 0.3]); % 95th        
        
        % formatting
        set(gca, 'xTick', 0:20:100)
        xlabel('% of periods in which value is exceeded');
        ylabel('total flow (ML) over period of specified length'); 
        
        switch iDur
            case 1, title('(a) Annual values'); legend('historic data', 'stoch. 5th & 95th', 'stoch. median')
            case 2, title('(b) 5 year runs, overlapping'); 
        end
        
        clearvars -except hist stoch ReqdDurations NumDur NumReplicates overlapping stoch_rearranged
    end
    
    % print the figure
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 22 12]); 
    print('-painters','-depsc', 'out\figures\Figure9'); saveas(gcf, 'out\figures\Figure9', 'png');   
    close
    
end

function Figure11(StochPertData, info, HistRecordLength)
    
    figure
    
    for iRow = 1:6
        
        spd = StochPertData; l = true;
        switch iRow
            case 1, lbl = 'a'; data.LessRisky = spd.deltaP_up;          data.MoreRisky = spd.deltaP_down; 
            case 2, lbl = 'b'; data.LessRisky = []; l = false;          data.MoreRisky = spd.deltaT_up;
            case 3, lbl = 'c'; data.LessRisky = spd.deltaLowFreqP_down; data.MoreRisky = spd.deltaLowFreqP_up;  
            case 4, lbl = 'd'; data.LessRisky = spd.deltaSeas_down;     data.MoreRisky = spd.deltaSeas_up; 
            case 5, lbl = 'e'; data.LessRisky = spd.deltaRRrel_down;    data.MoreRisky = spd.deltaRRrel_up;
            case 6, lbl = 'f'; data.LessRisky = []; l = false;          data.MoreRisky = spd.combination; 
        end
        
        %% column (i): seasonality
        
        % create an additional column in the tables containing the sum of the four unreg tribs - columns 4:7
        if iRow == 1, StochPertData.BaseCase.TS_Q_ReachInflows = AddCol_SumUnreg(StochPertData.BaseCase.TS_Q_ReachInflows); end
        if l, data.LessRisky.TS_Q_ReachInflows                 = AddCol_SumUnreg(data.LessRisky.TS_Q_ReachInflows); end
        data.MoreRisky.TS_Q_ReachInflows                       = AddCol_SumUnreg(data.MoreRisky.TS_Q_ReachInflows); 
        
        % convert the tables to monthly %
        data_i.StochBase.TS_Q       = ConvertToMonthlyPercent(StochPertData.BaseCase.TS_Q_ReachInflows);
        if l, data_i.LessRisky.TS_Q = ConvertToMonthlyPercent(data.LessRisky.TS_Q_ReachInflows); end
        data_i.MoreRisky.TS_Q       = ConvertToMonthlyPercent(data.MoreRisky.TS_Q_ReachInflows);
        
        % data formatting
        NumCols = 0; bp_data = [];
        for iMonth = 1:12

            % develop index identifying data for this month
            ind_base            = data_i.StochBase.TS_Q.Month == iMonth; NumMon_base   = sum(ind_base); 
            if l, ind_pert_less = data_i.LessRisky.TS_Q.Month == iMonth; end
            ind_pert_more       = data_i.MoreRisky.TS_Q.Month == iMonth; NumMon_pert  = sum(ind_pert_more); 
            
            % get data for this month - note, we are summing the unreg
            % inflows and thus summing columns 4:7
            ThisMon_base            = data_i.StochBase.TS_Q.SumUnreg(ind_base); 
            if l, ThisMon_pert_less = data_i.LessRisky.TS_Q.SumUnreg(ind_pert_less); end
            ThisMon_pert_more       = data_i.MoreRisky.TS_Q.SumUnreg(ind_pert_more);  
            
            % create nan-filled columns for data to go into
            bp_data (1:NumMon_base, (NumCols+1):(NumCols+4)) = nan; 

            % add data to boxplot table
            bp_data       (1:NumMon_base, NumCols+1) = ThisMon_base; 
            if l, bp_data (1:NumMon_pert, NumCols+2) = ThisMon_pert_less; end
            bp_data       (1:NumMon_pert, NumCols+3) = ThisMon_pert_more; 
            bp_ColourGroup    ((NumCols+1):(NumCols+4)) = [2 3 4 1];

            NumCols = NumCols + 4; 

        end
        
        % plotting
        p = 1+4*(iRow-1);
        subplot(6, 4, p:p+1); hold on;
        bp_labels = {'J', '', '', '', 'F', '', '', '', 'M', '', '', '', 'A', '', '', '', 'M', '', '', '', 'J', '', '', '', 'J', '', '', '', 'A', '', '', '', 'S', '', '', '', 'O', '', '', '', 'N', '', '', '', 'D', '', '', ''}; 
        colours_array = [  0   0   0; ... % blank series - colour irrelevant
                         0.5 0.5 0.5; ... % baseline - grey
                         0.3 0.3 0.8; ... %  risk - blue
                         0.8 0.3 0.3];    %  risk - red
        bp = boxplot(bp_data, 'PlotStyle', 'compact', 'ColorGroup', bp_ColourGroup, 'Labels', bp_labels, 'Colors', colours_array, 'Symbol', '');
        
        % fiddle with whisker length - make them 5th and 95th
        WhiskPct_lwr = 5; WhiskPct_upp = 95; 
        lines = get(gca, 'children');
        uw = findobj(lines, 'tag', 'Whisker'); numBx = size(uw, 1); assert(numBx == size(bp_data, 2));
        for i = 1:numBx, uw(i).YData = [prctile(bp_data(:, numBx-i+1), WhiskPct_lwr) prctile(bp_data(:, numBx-i+1), WhiskPct_upp)]; end
        
        % ylabel('Q in month (%)')
        ylabel({'these'; 'words'; 'are intended'; 'to take'; 'up'; 'space'}) % this is merely intended to shift the plot over a little; it will be deleted in post processing
        % title(['(' lbl '-i)']); 
        ylim([0 40])

        %% columns (ii) annual exceedence curves; and (iii) 5-year exceedence curves
        
        clearvars -except StochPertData info HistRecordLength lbl data iRow l
        
        % aggregate data to annual
        data.MoreRisky = AggregateToAnnual_stoch(data.MoreRisky, info.pars);
        if l, data.LessRisky = AggregateToAnnual_stoch(data.LessRisky, info.pars); end
        if iRow == 1, data.BaseCase = AggregateToAnnual_stoch(StochPertData.BaseCase, info.pars); end

        % get flow timeseries data to be plotted
        if l, iMax = 3; else, iMax = 2; end
        for i = 1:iMax
            names = {'BaseCase','MoreRisky','LessRisky'};
            flow_ts.(names{i}) = ...
               data.(names{i}).TS_annual_Q_ReachInflows.Reach1 ...
             + data.(names{i}).TS_annual_Q_ReachInflows.Reach2 ...
             + data.(names{i}).TS_annual_Q_ReachInflows.Reach3 ...
             + data.(names{i}).TS_annual_Q_ReachInflows.Reach4; 
        end
        
        % collate the data for the various options on show
        if l, MaxOption = 6; else, MaxOption = 4; end
        for iOption = 1:MaxOption
            
            clearvars ThisData stoch_rearranged StochRuns StochPctiles
            
            % define options
            switch iOption
                case 1, ThisDur = 1; ThisData = flow_ts.BaseCase;  ThisColour = [0.5 0.5 0.5]; ymax = 6; 
                case 2, ThisDur = 5; ThisData = flow_ts.BaseCase;  ThisColour = [0.5 0.5 0.5]; ymax = 20; 
                case 3, ThisDur = 1; ThisData = flow_ts.MoreRisky; ThisColour = [0.8 0.3 0.3]; ymax = 6;  
                case 4, ThisDur = 5; ThisData = flow_ts.MoreRisky; ThisColour = [0.8 0.3 0.3]; ymax = 20; 
                case 5, ThisDur = 1; ThisData = flow_ts.LessRisky; ThisColour = [0.3 0.3 0.8]; ymax = 6; 
                case 6, ThisDur = 5; ThisData = flow_ts.LessRisky; ThisColour = [0.3 0.3 0.8]; ymax = 20; 
            end
            
            % collate data
            stoch_rearranged = RearrangeData(HistRecordLength, ThisData); 
            NumReplicates = size(stoch_rearranged, 2); 

            % collate frequency curves
            overlapping = true; % periods can overlap
            sorted = true;
            for iRep = 1:NumReplicates
                StochRuns(:, iRep) = CollateFreqCurve(stoch_rearranged(:, iRep), ThisDur, overlapping, sorted); 
            end

            % extract various percentiles out of the stochastic runs
            percentiles = [1 5 10 20 50 70 90 95 99]; StochPctiles = [];
            for iPct = 1:size(percentiles, 2)
                NumOrd = size(StochRuns, 1); 
                for iOrd = 1:NumOrd
                    StochPctiles(iOrd, iPct) = prctile(StochRuns(iOrd, :), percentiles(iPct)); 
                end
            end            

            % choose which subplot for this option
            if ThisDur == 1, subplot(6, 4, 3+4*(iRow-1)); end
            if ThisDur == 5, subplot(6, 4, 4+4*(iRow-1)); end
            hold on
            
            % get x values
            clearvars x y
            n = size(StochPctiles, 1);
            x = 100 * ((1:n) / n); x = fliplr(x);
            
            % plot
            pos = 2; y(1:n) = 0.000001*StochPctiles(1:n, pos); plot(x, y, '--', 'Color', [ThisColour 0.6]); % 5th - note the last number is the transparency
            pos = 5; y(1:n) = 0.000001*StochPctiles(1:n, pos); plot(x, y, '-',  'Color', [ThisColour 1.0]); % median
            pos = 8; y(1:n) = 0.000001*StochPctiles(1:n, pos); plot(x, y, '--', 'Color', [ThisColour 0.6]); % 95th
            ylim([0 ymax]);
            
        end    
    end
   
    % print the figure
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 23 30]); 
    print('-painters','-depsc', 'out\figures\Figure11_raw'); saveas(gcf, 'out\figures\Figure11_raw', 'png');   
    close    
    
end

function Figure12()

    % Create a plot equivalent to Figure 12 from the paper

    % Note, the intention here is to demonstrate, rather than exactly
    % replicate, Figure 12 from the paper.  In line with this, we use dummy
    % data here rather than providing the actual data plotted in the paper.
    
    % load GCM data
    ClimProj = GetClimProj(); 
    
    for iFig = 1:2 
        
        clearvars -except iFig ClimProj
        close all;
    
        % define details and settings for each figure
        switch iFig
            case 1
                FileName = 'outTable_All_Alloc_mean.mat'; 
                settings.focus = 'HRWS'; % focus: high reliability water shares
                settings.zlimits = [0 1];
                settings.cm = flipud(parula); % colour map
                Outname = 'Figure12a'; 
            case 2
                FileName = 'outTable_All_stressConsec_mean.mat'; 
                settings.focus = 'LBF'; % focus: large bodied fish
                settings.zlimits = [0 20];
                settings.cm = copper;
                Outname = 'Figure12b'; 
        end
        settings.labels = true; % turn off for clean plots
        
        % load data
        data_all = GetData(FileName, settings.focus);
        
        % define axes for subplots
        PlotAxes (1, 1:4) = {'PvsT', 'PvsL', 'PvsS', 'PvsR'};
        PlotAxes (2, 1:4) = { 'n/a', 'TvsL', 'TvsS', 'TvsR'};
        PlotAxes (3, 1:4) = { 'n/a',  'n/a', 'LvsS', 'LvsR'};
        PlotAxes (4, 1:4) = { 'n/a',  'n/a',  'n/a', 'SvsR'};
        
        % loop to create each subplot
        figure
        for iRow = 1:4
            for iCol = 1:4
                % plot subplot
                CreateSubplot(iRow, iCol, PlotAxes{iRow, iCol}, data_all, settings, ClimProj); 
            end
        end
        
        % save plot
        set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 18 18]); 
        print('-painters','-depsc', ['out\figures\' Outname]); saveas(gcf, ['out\figures\' Outname], 'png');   
        close
    end    

end

function SuppMattPlots(HistoricData, StochPertData, info)
    
    % produce plots for Supplementary Material section S8 ("Plots for evaluation of quality of stochastic and perturbed data")

    % Figure S8.1: scatter plots of *monthly* T vs P, by season, for baseline stochastic data 
    mode = 'baseline';
    FigureS8_1(HistoricData, StochPertData, mode)

    % Figure S8.2: spatial dependencies between regions for annual P data
    FigureS8_2(HistoricData, StochPertData)
    
    % Figure S8.3: spatial dependencies between regions for monthly P data
    FigureS8_3(HistoricData, StochPertData)

    % Figure S8.4: spatial dependencies between regions for 5-year runs of P 
    FigureS8_4(HistoricData, StochPertData)

    % Figure S8.5: spatial dependencies between regions for annual Q data 
    FigureS8_5(HistoricData, StochPertData)

    % Figure S8.6: spatial dependencies between regions for monthly Q data 
    FigureS8_6(HistoricData, StochPertData)

    % Figure S8.7: spatial dependencies between regions for 5-year runs of Q
    FigureS8_7(HistoricData, StochPertData) 

    % Figure S8.8: spatial dependencies between regions for annual Q data, looking at perturbed stoch data 
    FigureS8_8(HistoricData, StochPertData)

    % Figure S8.9: spatial dependencies between regions for monthly Q data, looking at perturbed stoch data 
    FigureS8_9(HistoricData, StochPertData)

    % Figure S8.9: spatial dependencies between regions for monthly Q data, looking at perturbed stoch data 
    FigureS8_9(HistoricData, StochPertData)

    % Figure S8.10: scatter plots of *monthly* T vs P, by season, for perturbed stochastic data 
    mode = 'perturbed';
    FigureS8_1(HistoricData, StochPertData, mode)

    % Figure S8.11: plot of shifting of the P vs T relationship
    FigureS8_11(HistoricData, StochPertData)    

    % Figure S8.12 and S8.13: (i) Qave vs 1 in 20 driest 5 year runs; (ii) Qave vs seasonality
    HistRecordLength = size(HistoricData.flow_annual.ReachInflows_ML, 1);
    FigureS8_12_13(StochPertData, info, HistRecordLength)

    % Figure S8.14: QQ plots for historic and stochastic (incl. perturbed)
    FigureS8_14(HistoricData, StochPertData)

end

function StochTSdata = AggregateToAnnual_stoch(StochTSdata, pars)

    % aggregate the given stochastic monthly data to annual. 
    
    % There is an existing function to do this for the historic data, but
    % we can't use it as is because the names are different.  As a
    % workaround, this function renames as required, runs the existing
    % function, and then renames back at the end. 
    
    data_in.precip_monthly                           = StochTSdata.TS_P; 
    data_in.T_monthly                                = StochTSdata.TS_T; 
    data_in.PET_monthly                              = StochTSdata.TS_PET; 
    data_in.flow_monthly.ReachInflows_ML             = StochTSdata.TS_Q_ReachInflows; 
    data_in.flow_monthly.RepresentativeCatchments_mm = StochTSdata.TS_Q_ReprCatchments; 
    
    pars.checking = false; 
    data_out = AggregateToAnnual(data_in, pars);
    
    StochTSdata.TS_annual_P                = data_out.precip_annual;
    StochTSdata.TS_annual_T                = data_out.T_annual;
    StochTSdata.TS_annual_PET              = data_out.PET_annual;
    StochTSdata.TS_annual_Q_ReachInflows   = data_out.flow_annual.ReachInflows_ML;
    StochTSdata.TS_annual_Q_ReprCatchments = data_out.flow_annual.RepresentativeCatchments_mm;    

end

function UpdatedTable = ConvertToMonthlyPercent(OldTable)
    a = table2array(OldTable);
    AllYears = unique(a(:, 1)); 
    NumYears = size(AllYears, 1); 
    for iYear = 1:NumYears
        ThisYear = AllYears(iYear);     % what year are we working with
        ind = a(:, 1) == ThisYear;     % which rows
        FlowTot_allSubareas = sum(a(ind, 3:end)); 
        b = a(ind, 3:end);
        b = b ./ FlowTot_allSubareas;
        a(ind, 3:end) = b;
    end
    a(:, 3:end) = a(:, 3:end)*100; % proportion to percent
    UpdatedTable = array2table(a, 'VariableNames', OldTable.Properties.VariableNames); 
end

function StochRearranged = RearrangeData(HistRecordLength, StochOrig)
        
    % find the number of replicates and pre-allocate/pre-populate arrays for speed 
    length_stoch = size(StochOrig, 1); 
    NumReplicates = floor(length_stoch/HistRecordLength); 

    % do the rearrangement
    endpos = 0; 
    for iReplicate = 1:NumReplicates
        startpos = endpos + 1; 
        endpos = startpos + HistRecordLength - 1; 
        StochRearranged(1:HistRecordLength, iReplicate) = StochOrig(startpos:endpos, :); 
    end

end

function runs = CollateFreqCurve(ts, runlength, overlapping, sorted)

    % take a timeseries and:
    % - if 'sorted' = true:  create a frequency curve for multiannual runs of length runlength
    % - if 'sorted' = false: same again, but ordered according to their original sequence in time 
    % Runs can be either overlapping or not overlapping.
    
    pos = size(ts, 1); % start at the end and go backwards (for historic data, later periods 
                       % are considered more certain so we wouldn't want them ignored due to 
                       % being part of the unused remainder)                    

    % collate the data
    NumRuns = 0; 
    while pos-runlength > 0
        NumRuns = NumRuns + 1; 
        runs(NumRuns) = sum(ts((pos-runlength+1):pos));
        if overlapping
            pos = pos - 1; 
        else
            pos = pos - runlength; 
        end
    end

    % sort the data
    if sorted
        runs = sort(runs); 
    end
end

function outtable = AddCol_SumUnreg(intable)
    
    a = table2array(intable); 
    d = sum(a(:, 4:7), 2); % for this example, the unreg reaches are in positions 4 through 7 in intable
    outtable = array2table([a d], 'VariableNames', [intable.Properties.VariableNames'; 'SumUnreg']); 
    
end

function ClimProj = GetClimProj()
    
    % note: you need to go source the climate projections data, populate a file with it, and
    % then load this data, like this: 
    % ClimProj = load(<your file path>);

    % instead, so that you can run Figure 12, we produce some randomly generated dummy data:
    ave =  -3; std = 6.1; deltaP    = ave + std*randn(38, 1);
    ave =   1.49; std =  0.33; deltaT    = ave + std*randn(38, 1);
    ave = 0.0067; std = 0.018; deltaSeas = ave + std*randn(38, 1);
    DummyData = table(deltaP, deltaT, deltaSeas); % note, each row of this table is meant to represent a single GCM

    ClimProj.RCP85.Year_2040 = DummyData;

end

function data_all = GetData(FileName, focus)
    
    % note: you need to run the stochastic perturbed data through your
    % system model and calculate the system performance.  Then save this
    % information to file and load it here, like this: 
    % data_all = load([<your folder path> FileName]);    

    % instead, so that you can run Figure 12, we produce some randomly generated dummy data:
    deltaP_space           = [-0.40 -0.35 -0.30 -0.25 -0.20 -0.15 -0.10 -0.05 0 .05 .10 .15]; % rainfall proportional change
    deltaT_space           = [0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0];                             % additional degrees of warming
    deltaLowFreqP_space    = [-0.03 -0.015 0 +0.015 +0.03 +0.045 +0.06 +0.075];               % changes n Hurst Coefficient
    deltaSeasonality_space = [-0.06 -0.03 0.00 +0.03 +0.06 +0.09 +0.12 +0.15];                % changes in seasonality
    deltaRRrship_space     = [-50 -43.75 -37.5 -31.25 -25 -18.75 -12.5 -6.25 0 6.25 12.5];    % shift in rainfall runoff relationship
    
    % generate a table of nan's to hold the proxy data
    NumRows = size(deltaP_space, 2) * size(deltaT_space, 2) * size(deltaLowFreqP_space, 2) * size(deltaSeasonality_space, 2) * size(deltaRRrship_space, 2); 
    nan_array = nan(NumRows, 6);
    dummy_data_all = array2table(nan_array);
    dummy_data_all.Properties.VariableNames = {'deltaP', 'deltaT', 'deltaLowFreqP', 'deltaSeas', 'deltaRRrship', focus};

    % generate data
    i = 0; 
    for deltaP = deltaP_space
        a = (deltaP - min(deltaP_space))/range(deltaP_space);
        for deltaT = deltaT_space
            b = (deltaT - min(deltaT_space))/range(deltaT_space);
            for deltaLowFreqP = deltaLowFreqP_space
                c = (deltaLowFreqP - min(deltaLowFreqP_space))/range(deltaLowFreqP_space);
                for deltaSeasonality = deltaSeasonality_space    
                    d = (deltaSeasonality - min(deltaSeasonality_space))/range(deltaSeasonality_space);
                    for deltaRRrelship = deltaRRrship_space
                        e = (deltaRRrelship - min(deltaRRrship_space))/range(deltaRRrship_space);
                        
                        % calculate proxy values for performance
                        ThisPerf_dummy = 0.4 + a*0.4 - b*0.1 - c*0.05 - d*0.05 - e*0.25 + 0.1*rand();

                        % one more thing: in the paper, deltaPs are shown as percentages
                        deltaP_alt = 100*deltaP; 

                        % add values to table
                        i = i + 1; 
                        dummy_data_all{i, 1:6} = [deltaP_alt deltaT deltaLowFreqP deltaSeasonality deltaRRrelship ThisPerf_dummy];
                        
                    end
                end
            end
        end
    end
    
    % just a reminder that this is dummy data!
    data_all = dummy_data_all; 

end

function CreateSubplot(row, col, pa, data_all, settings, ClimProj)

    if ~strcmp(pa, 'n/a')
        
        % get x and y data
        for iDim = 1:2
            if iDim == 1, type = pa(1); else, type = pa(4); end
            switch type
                case 'P', xy.lab = [char(hex2dec('0394')) 'P'];     xy.data = data_all.("deltaP");         xy.gap = 0.050; xy.GCM = 1-ClimProj.RCP85.Year_2040{:, 1};
                case 'T', xy.lab = [char(hex2dec('0394')) 'T'];     xy.data = data_all.("deltaT");         xy.gap = 0.500; xy.GCM =   ClimProj.RCP85.Year_2040{:, 2};
                case 'L', xy.lab = [char(hex2dec('0394')) 'P-LF'];  xy.data = data_all.("deltaLowFreqP");  xy.gap = 0.015; xy.GCM =   nan(1, 38);
                case 'S', xy.lab = [char(hex2dec('0394')) 'S'];     xy.data = data_all.("deltaSeas");      xy.gap = 0.030; xy.GCM =   ClimProj.RCP85.Year_2040{:, 3};
                case 'R', xy.lab = [char(hex2dec('0394')) 'RRrel']; xy.data = data_all.("deltaRRrship");   xy.gap = 0.250; xy.GCM =   nan(1, 38);
            end           
            xy.uniq = unique(xy.data);
            xy.inc  = size(xy.uniq, 1);
            if iDim == 1, y = xy; else, x = xy; end
        end
        
        % get z values
        [x, y, z] = GetZvals(x, y, data_all, settings.focus);
        
        % create subplot for this plot
        subplot(4, 4, col + 4*(row - 1)); hold on
        colormap(settings.cm)
        p = pcolor(x.uniq, y.uniq, z');
        
        % add GCM scatter plot
        scatter(x.GCM, y.GCM, 15, 'w')

        % format plot
        xlim([min(x.uniq) max(x.uniq)]); ylim([min(y.uniq) max(y.uniq)]);
        set(p, 'EdgeColor', 'none');
        % caxis(settings.zlimits) % turn this on to control your colour bar limits
        if ~settings.labels, set(gca, 'visible', 'off'); end
        if row == 1, set(gca, 'Ydir', 'reverse'); end

        if settings.labels
            xlabel(x.lab); ylabel(y.lab);
            if(row*col == 1), colorbar; end
        else
            set(gca,'Yticklabel', [], 'Xticklabel', []) 
        end
    end

end

function [x, y, z] = GetZvals(x, y, data_all, focus)
    z = nan(x.inc, y.inc);
    for ix = 1:x.inc
        for iy = 1:y.inc

            % identify which rows of the input table are relevant to this combination
            ind_x = x.data == x.uniq(ix);
            ind_y = y.data == y.uniq(iy);
            ind = (ind_x + ind_y) == 2;

            % get z distribution for this combination of x and y
            z_dist = data_all.(focus)(ind); 

            % extract relevant statistic to plot
            z(ix, iy) = median(z_dist);

        end
    end

    % bewilderingly, and as Mathworks admit in their online help
    % for function "pcolor", 'none of the values in the last row or
    % column are represented in the plot'.  Go figure.  Therefore,
    % repeat the last row and column of the data.  
    x.uniq(x.inc+1) = x.uniq(x.inc) + x.gap; 
    z(x.inc+1, :) = nan;
    y.uniq(y.inc+1) = y.uniq(y.inc) + y.gap; 
    z(:, y.inc+1) = nan;
end

function FigureS8_1(HistoricData, StochPertData, mode)
    
    figure; reg = 'B'; % focus on region B
    for iSeason = 1:4

        switch iSeason
            case 1, s.name = 'Summer (DJF)'; m_ind = [12  1  2]; leg = {'Dec', 'Jan', 'Feb'};
            case 2, s.name = 'Autumn (MAM)'; m_ind = [ 3  4  5]; leg = {'Mar', 'Apr', 'May'};
            case 3, s.name = 'Winter (JJA)'; m_ind = [ 6  7  8]; leg = {'Jun', 'Jul', 'Aug'};
            case 4, s.name = 'Spring (SON)'; m_ind = [ 9 10 11]; leg = {'Sep', 'Oct', 'Nov'};
        end
        
        
        for iMonth = 1:3
            % get historic data for this season and month
            m = HistoricData.precip_monthly.Month; 
            ind = find(m==m_ind(iMonth)); sz_hist = size(ind, 1);
            histP  (:, iMonth) = HistoricData.precip_monthly.(reg)(ind); 
            histPET(:, iMonth) = HistoricData.PET_monthly.(reg)(ind);
            
            switch mode
                case 'baseline'
                    % get baseline stochastic data for this season
                    m = StochPertData.BaseCase.TS_P.Month;
                    ind = find(m==m_ind(iMonth)); 
                    ind = ind(1:(10*sz_hist), 1); 
                    stochP  (:, iMonth)   = StochPertData.BaseCase.TS_P.(reg)(ind);
                    stochPET(:, iMonth) = StochPertData.BaseCase.TS_PET.(reg)(ind); 

                    filepath = 'out\figures\FigureS8_1_raw';
                case 'perturbed'
                    % get perturbed stochastic data for this season
                    m = StochPertData.combination2.TS_P.Month;
                    ind = find(m==m_ind(iMonth)); 
                    ind = ind(1:(10*sz_hist), 1); 
                    stochP  (:, iMonth)   = StochPertData.combination2.TS_P.(reg)(ind);
                    stochPET(:, iMonth) = StochPertData.combination2.TS_PET.(reg)(ind);
                    
                    filepath = 'out\figures\FigureS8_10_raw';
            end
        end

        mkr = {'s', 'd', 'o'};
        for iPlot = 1:3
            subplot(4, 3, 3*(iSeason-1)+iPlot); hold on;
            
            % plot stochastic data in columns 2 and 3
            if iPlot ~= 1, for iMonth = 1:3, scatter(stochP(:, iMonth), stochPET(:, iMonth), 25, [0.85 0.3250 0.0980], 'marker', mkr{iMonth}); end; end
            
            % plot historic data in columns 1 and 3
            if iPlot ~= 2, for iMonth = 1:3, scatter(histP (:, iMonth),  histPET(:, iMonth), 25, [0.00 0.4470 0.7410], 'marker', mkr{iMonth}); end; end

            % plot legend for columns 1 & 2
            % if iPlot ~= 3, legend(leg); end

            % other formatting
            xlim([0 400]); ylim([0 200])
            
        end

    end

    % print the figure
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 17 20]); 
    print('-painters','-depsc', filepath); saveas(gcf, filepath, 'png');   
    close         
    
end

function FigureS8_2(HistoricData, StochPertData)
    
    % plot of spatial dependencies of annual precipitation

    % prepare inputs to bespoke function
    cols = [3; 6; 8; 9]; % corresponds to regions A, D, F and G
    data_hist  = HistoricData.precip_monthly{:, cols};
    data_stoch = StochPertData.BaseCase.TS_P{:, cols}; 
    xnam = {'Region A', 'Region D', 'Region F', 'Region G'};
    fpath = 'out\figures\FigureS8_2';
    AggMonToAnn = true; reversed = false; 

    % call bespoke function to create and save the figure
    kf_gplotmatrix(data_hist, data_stoch, xnam, fpath, AggMonToAnn, reversed)

end

function FigureS8_3(HistoricData, StochPertData)

    % plot of spatial dependencies of monthly precipitation
    
    % prepare inputs to bespoke function
    cols = [3; 6; 8; 9]; % corresponds to regions A, D, F and G
    data_hist  = HistoricData.precip_monthly{:, cols};
    data_stoch = StochPertData.BaseCase.TS_P{:, cols}; 
    xnam = {'Region A', 'Region D', 'Region F', 'Region G'};
    fpath = 'out\figures\FigureS8_3';
    AggMonToAnn = false; reversed = false; 

    % call bespoke function to create and save the figure
    kf_gplotmatrix(data_hist, data_stoch, xnam, fpath, AggMonToAnn, reversed)

end

function FigureS8_4(HistoricData, StochPertData)

    % plot of spatial dependencies of 5-year runs of precipitation
    
    % specify data sources
    cols = [3; 6; 8; 9]; % corresponds to regions A, D, F and G
    data_hist  = HistoricData.precip_monthly{:, cols};
    data_stoch = StochPertData.BaseCase.TS_P{:, cols};     
    
    % site-by-site processing
    runlength = 5; overlapping = true; sorted = false;
    for i = 1:4
        % for each site, aggregate to annual and get rolling sums
        h = sum(reshape(data_hist (:, i), 12, 107)); 
        s = sum(reshape(data_stoch(:, i), 12, 3000)); 
        runs_h(:, i) = CollateFreqCurve(h', runlength, overlapping, sorted);
        runs_s(:, i) = CollateFreqCurve(s', runlength, overlapping, sorted);
    end    

    % prepare remaining inputs to bespoke function
    xnam = {'Region A', 'Region D', 'Region F', 'Region G'};
    fpath = 'out\figures\FigureS8_4';
    AggMonToAnn = false; reversed = false; 

    % call bespoke function to create and save the figure
    kf_gplotmatrix(runs_h, runs_s, xnam, fpath, AggMonToAnn, reversed)

end

function FigureS8_5(HistoricData, StochPertData)
    
    % plot of spatial dependencies of annual streamflow

    % prepare inputs to bespoke function
    cols = 3:7; % corresponds to all inflows (Eildon plus reaches 1 - 4)
    data_hist  = 0.001 * HistoricData.flow_monthly.ReachInflows_ML{:, cols}; % note now GL not ML
    data_stoch = 0.001 * StochPertData.BaseCase.TS_Q_ReachInflows{:, cols}; % note now GL not ML
    xnam = {'Eildon inflow', 'Reach 1 inflow', 'Reach 2 inflow', 'Reach 3 inflow', 'Reach 4 inflow'};
    fpath = 'out\figures\FigureS8_5';
    AggMonToAnn = true; reversed = false; 

    % call bespoke function to create and save the figure
    kf_gplotmatrix(data_hist, data_stoch, xnam, fpath, AggMonToAnn, reversed)

end

function FigureS8_6(HistoricData, StochPertData)
    
    % plot of spatial dependencies of monthly streamflow

    % prepare inputs to bespoke function
    cols = 3:7; % corresponds to all inflows (Eildon plus reaches 1 - 4)
    data_hist  = 0.001 * HistoricData.flow_monthly.ReachInflows_ML{:, cols}; % note now GL not ML
    data_stoch = 0.001 * StochPertData.BaseCase.TS_Q_ReachInflows{:, cols};  % note now GL not ML
    xnam = {'Eildon inflow', 'Reach 1 inflow', 'Reach 2 inflow', 'Reach 3 inflow', 'Reach 4 inflow'};
    fpath = 'out\figures\FigureS8_6';
    AggMonToAnn = false; reversed = false; 

    % call bespoke function to create and save the figure
    kf_gplotmatrix(data_hist, data_stoch, xnam, fpath, AggMonToAnn, reversed)

end

function FigureS8_7(HistoricData, StochPertData)

    % plot of spatial dependencies of 5-year runs of precipitation
    
    % specify data sources
    cols = 3:7; % corresponds to all inflows (Eildon plus reaches 1 - 4)
    data_hist  = 0.001 * HistoricData.flow_monthly.ReachInflows_ML{:, cols}; % note now GL not ML
    data_stoch = 0.001 * StochPertData.BaseCase.TS_Q_ReachInflows{:, cols};  % note now GL not ML    
    
    % site-by-site processing
    runlength = 5; overlapping = true; sorted = false;
    for i = 1:5
        % for each site, aggregate to annual and get rolling sums
        h = sum(reshape(data_hist (:, i), 12, 127)); 
        s = sum(reshape(data_stoch(:, i), 12, 3000)); 
        runs_h(:, i) = CollateFreqCurve(h', runlength, overlapping, sorted);
        runs_s(:, i) = CollateFreqCurve(s', runlength, overlapping, sorted);
    end    

    % prepare remaining inputs to bespoke function
    xnam = {'Eildon inflow', 'Reach 1 inflow', 'Reach 2 inflow', 'Reach 3 inflow', 'Reach 4 inflow'};
    fpath = 'out\figures\FigureS8_7';
    AggMonToAnn = false; reversed = false; 

    % call bespoke function to create and save the figure
    kf_gplotmatrix(runs_h, runs_s, xnam, fpath, AggMonToAnn, reversed)

end

function FigureS8_8(HistoricData, StochPertData)
    
    % plot of spatial dependencies of annual streamflow

    % prepare inputs to bespoke function
    cols = 3:7; % corresponds to all inflows (Eildon plus reaches 1 - 4)
    data_hist  = 0.001 * HistoricData.flow_monthly.ReachInflows_ML{:, cols}; % note now GL not ML
    data_stoch = 0.001 * StochPertData.combination2.TS_Q_ReachInflows{:, cols}; % note now GL not ML
    xnam = {'Eildon inflow', 'Reach 1 inflow', 'Reach 2 inflow', 'Reach 3 inflow', 'Reach 4 inflow'};
    fpath = 'out\figures\FigureS8_8';
    AggMonToAnn = true; reversed = true; 

    % call bespoke function to create and save the figure
    kf_gplotmatrix(data_hist, data_stoch, xnam, fpath, AggMonToAnn, reversed)

end

function FigureS8_9(HistoricData, StochPertData)
    
    % plot of spatial dependencies of monthly streamflow

    % prepare inputs to bespoke function
    cols = 3:7; % corresponds to all inflows (Eildon plus reaches 1 - 4)
    data_hist  = 0.001 * HistoricData.flow_monthly.ReachInflows_ML{:, cols}; % note now GL not ML
    data_stoch = 0.001 * StochPertData.combination2.TS_Q_ReachInflows{:, cols};  % note now GL not ML
    xnam = {'Eildon inflow', 'Reach 1 inflow', 'Reach 2 inflow', 'Reach 3 inflow', 'Reach 4 inflow'};
    fpath = 'out\figures\FigureS8_9';
    AggMonToAnn = false; reversed = true; 

    % call bespoke function to create and save the figure
    kf_gplotmatrix(data_hist, data_stoch, xnam, fpath, AggMonToAnn, reversed)

end

function FigureS8_11(HistoricData, StochPertData)
    
    % plots of P vs T for each month
    
    % get stochastic data
    x_all = StochPertData.BaseCase.TS_P.D(1:(1000*12)); 
    y_all = StochPertData.BaseCase.TS_T.D(1:(1000*12)); 
    months = StochPertData.BaseCase.TS_P.Month(1:(1000*12)); 
    
    % get historic data
    x_hist = HistoricData.precip_monthly.D;
    y_hist = HistoricData.T_monthly.D;
    months_hist = HistoricData.T_monthly.Month; 

    % define titles for each subplot
    titles = {  '(a) Jan - hist. vs sto.','(b) Jan - stoch. pert.', ...
                '(c) Feb - hist. vs sto.','(d) Feb - stoch. pert.', ...
                '(e) Mar - hist. vs sto.','(f) Mar - stoch. pert.', ...
                '(g) Apr - hist. vs sto.','(h) Apr - stoch. pert.', ...
                '(i) May - hist. vs sto.','(j) May - stoch. pert.', ...
                '(k) Jun - hist. vs sto.','(l) Jun - stoch. pert.', ...
                '(m) Jul - hist. vs sto.','(n) Jul - stoch. pert.', ...
                '(o) Aug - hist. vs sto.','(p) Aug - stoch. pert.', ...
                '(q) Sep - hist. vs sto.','(r) Sep - stoch. pert.', ...
                '(s) Oct - hist. vs sto.','(t) Oct - stoch. pert.', ...
                '(u) Nov - hist. vs sto.','(v) Nov - stoch. pert.', ...
                '(w) Dec - hist. vs sto.','(x) Dec - stoch. pert.'};
    
    figure;
    for iMonth = 1:12

        % prep
        ind = months == iMonth;
        x = x_all(ind); y = y_all(ind);

        % plotting - historic data
        plotpos = iMonth*2-1;
        subplot(6, 4, plotpos); hold on
        scatter(-99, -99, 2, [0.00 0.75 0.88]); scatter(-99, -99, [0.85 0.3250 0.0980]); % for benefit of legend only
        scatter(x, y, 10, [0.85 0.3250 0.0980], 'o', 'MarkerEdgeAlpha', 0.1) % baseline
        x2 = x_hist(months_hist == iMonth); y2 = y_hist(months_hist == iMonth);
        scatter(x2, y2, 2, [0.00 0.75 0.88], 'o', 'filled') % historic
        xlim([0 500]); ylim([0 40])
        text(510, 40, titles{plotpos}, 'FontSize',8,'HorizontalAlignment','right');

        % features that appear only once, such as legend and axis labels
        if iMonth == 1,  legend('hist.', 'stoch.', 'Location', 'southeast'); end 
        if iMonth == 9,  text(-170, 40, 'temperature (°C)', 'FontSize',16, 'Rotation',90); end
        if iMonth == 12, text(+500, -17, 'precipitation (mm/month)', 'FontSize',16, 'HorizontalAlignment','right'); end

        % plotting - perturbed data
        plotpos = iMonth*2;
        subplot(6, 4, plotpos); hold on
        scatter(-99, -99, 'yellow'); scatter(-99, -99, 'green'); scatter(-99, -99, 'magenta'); % for benefit of legend only
        scatter(x, y+4, 10, 'yellow', 'o', 'MarkerEdgeAlpha', 0.1) % higher temperature
        scatter(x*1.15, y, 10, 'green', 'o', 'MarkerEdgeAlpha', 0.1) % more rainfall
        scatter(x*0.7, y, 10, 'magenta', 'o', 'MarkerEdgeAlpha', 0.1) % less rainfall
        xlim([0 500]); ylim([0 40])
        text(510, 40, titles{plotpos}, 'FontSize',8,'HorizontalAlignment','right');
        if iMonth == 1, legend('T +4°C', 'P +15%', 'P -30% ', 'Location', 'southeast'); end
    end

    % print the figure
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 20 24]); 
    print('-painters','-depsc', 'out\figures\FigureS8_11'); saveas(gcf, 'out\figures\FigureS8_11', 'png');   
    close      

end

function FigureS8_12_13(StochPertData, info, HistRecordLength)
    
    % start two separate figures
    Fig5yr  = figure; % for the 5-year runs
    FigSumm = figure; % for the summer fraction
    
    for iPlot = 1:6
        
        spd = StochPertData; l = true;
        switch iPlot
            case 1, lbl = '(a) axis 1: mean annual precip.';     data.LessRisky = spd.deltaP_up;          data.MoreRisky = spd.deltaP_down; 
            case 2, lbl = '(b) axis 2: mean temperature';              data.LessRisky = []; l = false;          data.MoreRisky = spd.deltaT_up;
            case 3, lbl = '(c) axis 3: precipitation_{low frequency}'; data.LessRisky = spd.deltaLowFreqP_down; data.MoreRisky = spd.deltaLowFreqP_up;  
            case 4, lbl = '(d) axis 4: seasonality of precip.';  data.LessRisky = spd.deltaSeas_down;     data.MoreRisky = spd.deltaSeas_up; 
            case 5, lbl = '(e) axis 5: rainfall runoff rel''ship';  data.LessRisky = spd.deltaRRrel_down;    data.MoreRisky = spd.deltaRRrel_up;
            case 6, lbl = '(f) combination of small perturbations';    data.LessRisky = []; l = false;          data.MoreRisky = spd.combination; 
        end
        data.BaseCase = StochPertData.BaseCase;
        
        % for each of 23 sequences of 126 years: 

        % get average flows and 95th percentile 5-year runs
        [meanQ, FiveYr5thPct] = GetQannual_stats(data, l, HistRecordLength, info);        
        
        % get fraction of Q in summer
        frac_Q_summer = GetFracQ_summer(data, l, HistRecordLength);
        
        % Plot
        for iFig = 1:2
            
            % define options
            switch iFig
                case 1, ax = subplot(3, 2, iPlot, 'Parent', Fig5yr);  x = meanQ; y = FiveYr5thPct; m = 0.001; yl = [0 6000]; al1 = -1200; al2 = 6000; al3 = '5-year flow volume_{5th percentile}'; 
                case 2, ax = subplot(3, 2, iPlot, 'Parent', FigSumm); x = meanQ; y = frac_Q_summer; m = 1; yl = [0 0.4]; al1 = -0.1; al2 = 0.4; al3 = 'fraction of flow in warm months (NDJFMA)'; 
            end
            hold(ax, 'on')

            % plot
            scatter      (ax, 0.001*x.base, m*y.base, 25, [0.5000 0.5000 0.5000])   % baseline
            if l, scatter(ax, 0.001*x.less, m*y.less, 25, [0.0000 0.4470 0.7410]); end % less risky
            scatter      (ax, 0.001*x.more, m*y.more, 25, [0.8500 0.3250 0.0980]);  % more risky
            
            % formatting
            set(ax, 'xlim', [0 2000], 'ylim', yl)
            text(ax, 70, yl(2)*0.93, lbl, 'FontSize', 8)
            if iPlot == 5
                text(ax, 1500, al1, 'mean annual flow (GL)', 'FontSize', 14); 
                text(ax, -500, al2, al3, 'FontSize', 14, 'Rotation', 90); 
            end

            if iPlot == 1
                legend(ax, 'baseline', 'pert., risk \downarrow', 'pert., risk \uparrow', 'Location', 'southeast'); 
            end
        end
    end

    % print the figures
    set(Fig5yr, 'PaperUnits', 'centimeters'); set(Fig5yr, 'PaperPosition', [0 0 15 18]); 
    print(Fig5yr, '-painters','-depsc', 'out\figures\FigureS8_12'); saveas(Fig5yr, 'out\figures\FigureS8_12', 'png');   
    set(FigSumm, 'PaperUnits', 'centimeters'); set(FigSumm, 'PaperPosition', [0 0 15 18]); 
    print(FigSumm, '-painters','-depsc', 'out\figures\FigureS8_13'); saveas(FigSumm, 'out\figures\FigureS8_13', 'png');       
    
    close all;

end

function FigureS8_14(HistoricData, StochPertData)

    subarea = 'B'; lab = {'a','b','c','d'};
    
    % get annual historical data
    hist.annual = HistoricData.precip_annual.(subarea);

    % get annual stochastic data
    s = StochPertData; % note this is in monthly timestep
    d1 = s.BaseCase.          TS_P.(subarea); d2 = sum(reshape(d1, 12, (1/12)*size(d1, 1)))'; stoch.base.annual = d2; % middle step is aggregation to annual
    d1 = s.deltaLowFreqP_up.  TS_P.(subarea); d2 = sum(reshape(d1, 12, (1/12)*size(d1, 1)))'; stoch.up.annual   = d2;
    d1 = s.deltaLowFreqP_down.TS_P.(subarea); d2 = sum(reshape(d1, 12, (1/12)*size(d1, 1)))'; stoch.down.annual = d2;

    figure; iPlot = 0; 
    for iDur = [2 3 5 10]
        
        % put into runs
        runlength_yrs = iDur; overlapping = false; sorted = false;
        str = ['runs' num2str(runlength_yrs) 'yr'];
        d3 = hist.annual;       hist.      (str).data = CollateFreqCurve(d3, runlength_yrs, overlapping, sorted);
        d3 = stoch.base.annual; stoch.base.(str).data = CollateFreqCurve(d3, runlength_yrs, overlapping, sorted);
        d3 = stoch.up.  annual; stoch.up  .(str).data = CollateFreqCurve(d3, runlength_yrs, overlapping, sorted);
        d3 = stoch.down.annual; stoch.down.(str).data = CollateFreqCurve(d3, runlength_yrs, overlapping, sorted);

        % activate subplot
        iPlot = iPlot + 1; 
        subplot(2, 2, iPlot)
        
        % analyse using qq plot function (NOTE because hold is *not* on, the plot is getting wiped each time)
        h = qqplot(hist.(str).data, stoch.base.(str).data); stoch.base.(str).qqx = h(1).XData; stoch.base.(str).qqy = h(1).YData; 
        h = qqplot(hist.(str).data, stoch.up  .(str).data); stoch.up  .(str).qqx = h(1).XData; stoch.up  .(str).qqy = h(1).YData; 
        h = qqplot(hist.(str).data, stoch.down.(str).data); stoch.down.(str).qqx = h(1).XData; stoch.down.(str).qqy = h(1).YData; 

        % final plotting
        x = stoch.base.(str).qqx; y = stoch.base.(str).qqy; scatter(0.001*x, 0.001*y, 35, [0.20 0.2000 0.2000], 'o'); hold on
        x = stoch.up  .(str).qqx; y = stoch.up  .(str).qqy; scatter(0.001*x, 0.001*y, 35, [0.85 0.3250 0.0980], 'o', 'filled')
        x = stoch.down.(str).qqx; y = stoch.down.(str).qqy; scatter(0.001*x, 0.001*y, 35, [0.00 0.4470 0.7410], 'o', 'filled')
        yl = ylim; xl = xlim; % find current plot limits
        plot([0 1E10], [0 1E10], 'Color', [0.8 0.8 0.8], 'LineStyle','--', 'LineWidth', 2); % 1 to 1 line
        x = stoch.base.(str).qqx; y = stoch.base.(str).qqy; scatter(0.001*x, 0.001*y, 35, [0.20 0.2000 0.2000], 'o'); % repeat of base case to ensure it gets plotted on top
        
        % formatting
        ord_min = min([xl(1) yl(1)]); ord_max = max([xl(2) yl(2)]); % find plot limits
        ylim([ord_min ord_max]); xlim([ord_min ord_max]);
        x = ord_min + 0.05*(ord_max - ord_min);
        y = ord_max - 0.05*(ord_max - ord_min);
        text(x, y, ['(' lab{iPlot} ') ' num2str(runlength_yrs) ' year runs, non-overlapping'])
        if iPlot == 1, legend('baseline', 'Hurst +0.06', 'Hurst -0.015', 'Location','southeast'); end
        xlabel('quantile, historic (GL)'); ylabel('quantile, stochastic (GL)')
        grid on

    end

    % print the figure
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 18 18]); 
    print('-painters','-depsc', 'out\figures\FigureS8_14'); saveas(gcf, 'out\figures\FigureS8_14', 'png');   
    close       
    
end

function kf_gplotmatrix(data_hist, data_stoch, xnam, fpath, AggMonToAnn, reversed)
    
    NumSites = size(data_hist, 2);

    % if required, aggregate each location's annual timeseries from monthly to annual
    if AggMonToAnn
        
        % define number of years
        NumYears_hist  = (1/12)*size(data_hist,  1);
        NumYears_stoch = (1/12)*size(data_stoch, 1);        
        
        % aggregate to annual
        for i = 1:NumSites
            h(i, :) = sum(reshape(data_hist (:, i), 12, NumYears_hist)); 
            s(i, :) = sum(reshape(data_stoch(:, i), 12, NumYears_stoch)); 
        end
        data_hist  = h'; 
        data_stoch = s'; 
        n = NumYears_hist;
    else
        n = size(data_hist, 1);
    end

    % transpose and, in the case of stoch_final, shorten to 10x the historic length
    if size(data_stoch, 1)<10*n
        error('You need stochastic series at least 10 times as long as historic, to run this function.'); 
    end
    data_stoch = data_stoch(1:(10*n), :);

    % prepare labels
    [labels_1{1:n}] = deal('historic'); 
    [labels_2{1:(10*n)}] = deal('stoch.'); 

    % combine together
    if ~reversed
        data_all = [data_stoch; data_hist];
        labels_all = [labels_2'; labels_1'];
        clr = [0.85 0.3250 0.0980; 0.00 0.4470 0.7410]; siz = [2, 4.5];
    else
        data_all = [data_hist; data_stoch];
        labels_all = [labels_1'; labels_2'];
        clr = [0.00 0.4470 0.7410; 0.85 0.3250 0.0980]; siz = [4.5, 2];
    end


    % plot using gplotmatrix
    sym = 'oo'; doleg = 'on'; dispopt = []; 
    gplotmatrix(data_all, [], labels_all, clr, sym, siz, doleg, dispopt, xnam);
    

    % print the figure
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 17 16]); 
    warning('off','all')
    print('-painters','-depsc', fpath); saveas(gcf, fpath, 'png');   
    warning('on','all')
    close

end

function frac_Q_summer = GetFracQ_summer(AllData, l, HistRecordLength)
    
    if l, iMax = 3; else, iMax = 2; end % only do the "less risky" one if it exists
    
    for i = 1:iMax

        clearvars -except AllData l HistRecordLength i iMax frac_Q_summer

        switch i
            case 1, ThisData = AllData.BaseCase.TS_Q_ReachInflows; 
            case 2, ThisData = AllData.MoreRisky.TS_Q_ReachInflows;
            case 3, ThisData = AllData.LessRisky.TS_Q_ReachInflows;
        end

        % create an additional column in the tables containing the sum of the four unreg tribs - columns 4:7
        data = AddCol_SumUnreg(ThisData); 
        
        % data formatting
        for iMonth = 1:12
    
            % develop index identifying data for this month
            ind = data.Month == iMonth; NumMon = sum(ind);

            % get data for this month
            ThisMon = data.SumUnreg(ind);  
            
            % reshape into sets of 126 years (same as historic record length)
            n = floor(NumMon / HistRecordLength); m = n*HistRecordLength;
            ThisMon = reshape(ThisMon(1:m), HistRecordLength, n);
    
            % take averages for this month across each set of 126 years
            a(iMonth, 1:n) = sum(ThisMon);
    
        end
        
        % get summer fraction
        summer_ind = [1 2 3 4 11 12]; 
        this_frac_Q_summer = sum(a(summer_ind, :))./sum(a(:, :));
        
        % save info
        switch i
            case 1, frac_Q_summer.base = this_frac_Q_summer;
            case 2, frac_Q_summer.more = this_frac_Q_summer;
            case 3, frac_Q_summer.less = this_frac_Q_summer;
        end
    end
end

function [meanQ_all, FiveYr5thPct_all] = GetQannual_stats(AllData, l, HistRecordLength, info)
    
    if l, iMax = 3; else, iMax = 2; end % only do the "less risky" one if it exists
    
    for i = 1:iMax

        clearvars -except AllData l HistRecordLength info i iMax meanQ_all FiveYr5thPct_all

        switch i
            case 1, ThisData = AllData.BaseCase; 
            case 2, ThisData = AllData.MoreRisky;
            case 3, ThisData = AllData.LessRisky;
        end    

        % aggregate data to annual
        ThisData_agg = AggregateToAnnual_stoch(ThisData, info.pars);
        
        % get flow timeseries data to be plotted
        flow_ts = ThisData_agg.TS_annual_Q_ReachInflows.Reach1 ...
                + ThisData_agg.TS_annual_Q_ReachInflows.Reach2 ...
                + ThisData_agg.TS_annual_Q_ReachInflows.Reach3 ...
                + ThisData_agg.TS_annual_Q_ReachInflows.Reach4;
    
        % collate data
        stoch_rearranged = RearrangeData(HistRecordLength, flow_ts); 
        NumReplicates = size(stoch_rearranged, 2); 
    
        % collate frequency curves for 5 year overlapping runs
        overlapping = true; sorted = true; ThisDur = 5;
        for iRep = 1:NumReplicates
            StochRuns(:, iRep) = CollateFreqCurve(stoch_rearranged(:, iRep), ThisDur, overlapping, sorted); 
        end 
    
        % get meanQ
        meanQ = mean(stoch_rearranged); 
        
        % get 5th percentile 5-year flow (ie. lower than 95% of runs)
        for iRep = 1:NumReplicates
            FiveYr5thPct(iRep) = prctile(StochRuns(:, iRep), 5); 
        end

        % save info
        switch i
            case 1, meanQ_all.base = meanQ; FiveYr5thPct_all.base = FiveYr5thPct; 
            case 2, meanQ_all.more = meanQ; FiveYr5thPct_all.more = FiveYr5thPct; 
            case 3, meanQ_all.less = meanQ; FiveYr5thPct_all.less = FiveYr5thPct; 
        end        


    end
end
