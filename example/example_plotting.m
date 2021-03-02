
function example_plotting(StochPertData, HistoricData, info)

    % Created January 2021 by Keirnan Fowler, University of Melbourne, fowler.k@unimelb.edu.au

    % Integrated framework for rapid climate stress testing on a monthly timestep
    % by Keirnan Fowler, Natasha Ballis, Avril Horne, Andrew John, Rory Nathan and Murray Peel 

    % Licence: CC BY 3.0 - see https://creativecommons.org/licenses/by/3.0/au/

    % The purpose of this function is to produce the plots for the example
    % given in the paper.  Note: in the two more complex figures (Figure 6a 
    % and Figure 10) this code provides only the base data that is plotted 
    % and *not* extras such as figure labels, axis labels, and legends.  
    % This is to keep the code simple and also because the authors added 
    % these features in postprocessing using Inkscape (https://inkscape.org/).  
    
    if ~info.isOctave

%         % Figure 6a: separation of historic annual rainfall into low and high frequency components 
%         Figure7a(HistoricData, info);
%         disp('Figure 7a done.'); 
% 
%         % Figure 6b: combination of low frequency component across subareas
%         Figure7b(HistoricData, info); 
%         disp('Figure 7b done.'); 
% 
%         % Figure 7: seasonality of precipitation, temperature, PET and streamflow
%         Figure8(HistoricData, StochPertData);
%         disp('Figure 8 done.'); 

        % Figure 8: Match of stochastic baseline streamflow with historic for
        % annual and 5 year totals
        Figure9(HistoricData, StochPertData, info)
        disp('Figure 9 done.'); 

        % Figure 10: Evaluation of stochastic data
        HistRecordLength = size(HistoricData.flow_annual.ReachInflows_ML, 1);
        Figure11(StochPertData, info, HistRecordLength);
        disp('Figure 11 done.');
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
        HistRuns = CollateFreqCurve(hist, ThisDur, overlapping); 
        for iRep = 1:NumReplicates
            StochRuns(:, iRep) = CollateFreqCurve(stoch_rearranged(:, iRep), ThisDur, overlapping); 
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
            for iRep = 1:NumReplicates
                StochRuns(:, iRep) = CollateFreqCurve(stoch_rearranged(:, iRep), ThisDur, overlapping); 
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

function runs = CollateFreqCurve(ts, runlength, overlapping)

    % take a timeseries and create a frequency curve for multiannual runs
    % of length runlength, either overlapping or not overlapping.
    
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
    runs = sort(runs); 
end

function outtable = AddCol_SumUnreg(intable)
    
    a = table2array(intable); 
    d = sum(a(:, 4:7), 2); % for this example, the unreg reaches are in positions 4 through 7 in intable
    outtable = array2table([a d], 'VariableNames', [intable.Properties.VariableNames'; 'SumUnreg']); 
    
end