function [TS_P, extra_h] = PerturbP_Seasonality(TS_P_interim, deltaSeasonality, info)

    % Created 08/10/2019 by Keirnan Fowler, University of Melbourne
    
    % The purpose of this code is to perturb the seasonality of
    % precipitation through transferring some rainfall from the cold season
    % to the warm season (or vice-versa).  Prior to perturbation, seasonality
    % is based on historic seasonality via method of fragments.  Via 
    % this function, the precip. seasonality is changed so that more (less)
    % rainfall occurs in the warm season for positive (negative) values of
    % deltaSeasonality.  
    
    % for each subarea
    SubareaList = info.SubareaList'; TS_P = [TS_P_interim.Year TS_P_interim.Month]; CF_all = [];
    for subarea = SubareaList % equivalent of VBA "For Each subarea in SubareaList"
            
        % initialise vectors with required info
        ListOfYears       = TS_P_interim.Year;
        ListOfMonths      = TS_P_interim.Month;
        this_TS_P_interim = TS_P_interim.(subarea{:}); % unperturbed timeseries
        
        % initialise temporary vector for perturbed Precip
        P_pert = nan(size(this_TS_P_interim)); 

        % initial analysis to determine adjustment factors for each month
        [MonthlyAdjFactors, MonthlyArray_orig, MonthlyArray_target, extra_h_frag] = GetTargetDist(this_TS_P_interim, ListOfYears, ListOfMonths, deltaSeasonality); 
        extra_h.OrigMean.(subarea{:}) = extra_h_frag.OrigMean; 

        % loop over years
        for iYear = 1:info.pars.StochRepLen_yrs

            % pull out this year
            startmonth = (1+12*(iYear-1)); endmonth = 12*iYear;
            ind = startmonth:endmonth; 
            orig.precip = this_TS_P_interim(ind);
            orig.total  = sum(orig.precip); 

            % apply adjustment factors
            for iMonth = 1:12, pert.precip(iMonth) = orig.precip(iMonth) * MonthlyAdjFactors(iMonth); end

            % annual check to ensure the annual total is the same as what we started with 
            pert.total = sum(pert.precip); 
            CorrectionFactor(iYear) = orig.total / pert.total; 
            pert.precip = pert.precip * CorrectionFactor(iYear); 
            assert(abs(sum(pert.precip) - sum(orig.precip))<0.001);

            % add this back into the temporary array
            P_pert(ind) = pert.precip; 
        end

        % append to array
        TS_P = [TS_P P_pert]; 
        CF_all = [CF_all CorrectionFactor'];
        
        % % optional plotting
        % for iMonth = 1:12, AdoptedDist(iMonth, :) = TS_P_TempArray(ListOfMonths == iMonth); end
        % RegionName = region{1}; 
        % dummy = PlotSeasonalityPerturb1(MonthlyArray_orig, MonthlyArray_target, AdoptedDist, RegionName, deltaSeasonality); % monthly means and distributions
        % dummy = PlotSeasonalityPerturb2(this_TS_P_interim, TS_P_TempArray, RegionName, deltaSeasonality); % timeseries
    end
    
    % convert outputs to tables
    TS_P = array2table(TS_P, 'VariableNames', [{'Year'; 'Month'}; info.SubareaList]); 
    extra_h.CorrectionFactors = array2table(CF_all, 'VariableNames', info.SubareaList);
end

function [MonthlyAdjFactors, MonthlyArray_orig, MonthlyArray_target, extra_h_frag] = GetTargetDist(this_TS_P_interim, ListOfYears, ListOfMonths, deltaSeasonality)
    
    % get mean precipitation
    MeanAnnualP = sum(this_TS_P_interim)/(size(this_TS_P_interim, 1)/12);
    
    % determine total amount to be transferred from cooler months to warmer, long term average  
    transfer.total = deltaSeasonality * MeanAnnualP; % note, negative values mean a transfer in the opposite direction.

    % assign the seasonal pattern of transfer.  Notes: 1. these seasons are in line with the Climate Change in Australia website
    %                                                  2. the words 'gains' and 'loses' need to be switched if deltaSeasonality is negative 
    transfer.SeasPatt = [+0.2143  ... % January   - warm month, gains 3/14 of total transfer
                         +0.2143  ... % February  - warm month, gains 3/14 of total transfer
                         +0.2143  ... % March     - warm month, gains 3/14 of total transfer
                          0.0714  ... % April     - shoulder,   gains 1/14 of total transfer
                         -0.0714  ... % May       - shoulder,   loses 1/14 of total transfer
                         -0.2143  ... % June      - cold month, loses 3/14 of total transfer
                         -0.2143  ... % July      - cold month, loses 3/14 of total transfer
                         -0.2143  ... % August    - cold month, loses 3/14 of total transfer
                         -0.2143  ... % September - cold month, loses 3/14 of total transfer
                         -0.0714  ... % October   - shoulder,   loses 1/14 of total transfer
                         +0.0714  ... % November  - shoulder,   gains 1/14 of total transfer
                         +0.2143  ... % December  - warm month, gains 3/14 of total transfer
                         ];
    
    % % code to produce raw version of Figure S1:
    % plot(100*transfer.SeasPatt, '--x', 'MarkerSize', 12); xlim([0 13]); ylim([-30 30]); set(gca, 'XTick', 1:12);
    % set(gca, 'XTick', 1:12, 'XTickLabel', {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'});                     
    
    % assign the average amount (mm) of transfer for each month
    transfer.MonthlyAmounts = transfer.total * transfer.SeasPatt; 
    
    % for each month
    for iMonth = 1:12
        
        % get average P value for this month (in mm/month)
        ind = ListOfMonths == iMonth;
        OrigMean = mean(this_TS_P_interim(ind)); 
        
        % define mm amount as a proportion
        MonthlyAdjFactors(iMonth) = 1+(transfer.MonthlyAmounts(iMonth) / OrigMean);
        assert(min(MonthlyAdjFactors)>0, 'Error: you have pushed the seasonality perturbation too far: please revise with less extreme perturbation.'); 
        
        % define original distribution
        MonthlyArray_orig(iMonth, :) = this_TS_P_interim(ind);
        
        % define target distribution
        MonthlyArray_target(iMonth, :) = MonthlyArray_orig(iMonth, :) * MonthlyAdjFactors(iMonth); 
        
        % add variables to ancillary output extra_h_frag
        extra_h_frag.OrigMean(iMonth) = OrigMean; 
        
    end
    
end

function dummy = PlotSeasonalityPerturb1(OrigDist, TargetDist, AdoptedDist, RegionName, deltaSeasonality) 
        
    % plot 1 - means
    %%
    figure; 
    subplot(2, 1, 1); hold on;
    title(['(a) monthly means across all ' num2str(size(OrigDist, 2)) ' years in replicate'])
    plot(mean(OrigDist'), 'b-o');
    plot(mean(TargetDist'), '-x', 'color', [1 0.6 0.6], 'MarkerSize', 15)
    plot(mean(AdoptedDist'), 'r-o');
    xlim([0 12.99]); ylim([0 150])
    legend('original (synthetic unperturbed)', 'target', 'adopted')
    xlabels = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    set(gca, 'XTick', 1:12, 'XTickLabel', xlabels)
    ylabel('Monthly rainfall (mm/month)')
    grid on
    
    % plot 2 - boxplots, all on one figure
    subplot(2, 1, 2); hold on;
    title('(b) monthly distributions')
    
    % data prep for plot 2
    NumMonths = size(OrigDist, 2);
    BoxPlotData = NaN(NumMonths, 48); group = NaN(NumMonths, 48); pos = 0;
    for iMonth = 1:12
        pos = pos+1; BoxPlotData(:, pos) = OrigDist   (iMonth, :)'; group(:, pos) = 1*ones(NumMonths, 1)';
        pos = pos+1; BoxPlotData(:, pos) = TargetDist (iMonth, :)'; group(:, pos) = 1*ones(NumMonths, 1)';
        pos = pos+1; BoxPlotData(:, pos) = AdoptedDist(iMonth, :)'; group(:, pos) = 1*ones(NumMonths, 1)';
        pos = pos+1; % blank column
    end
    
    % plotting for plot 2
    b = boxplot(BoxPlotData, "PlotStyle", "compact"); 
    c = get(get(gca,'children'),'children');
    for i = 1:48
        set(c(i), 'String', ''); % labels - remove them
        set(c(i+1*48), 'MarkerEdgeColor', [0.9 0.9 0.9], 'MarkerSize', 2); % outliers
        switch mod(i, 4)
            case 0 % original (note we go backwards from here)
                set(c(i+2*48), 'MarkerEdgeColor', [0 0 1]); % median inner
                set(c(i+3*48), 'MarkerEdgeColor', [0 0 1]); % median outer
                set(c(i+4*48), 'color', [0 0 1]); % box
                set(c(i+5*48), 'color', [0 0 1]); % whisker               
            case 1 % no change - blank
            case 2 % adopted
                set(c(i+2*48), 'MarkerEdgeColor', [1 0 0]); % median inner
                set(c(i+3*48), 'MarkerEdgeColor', [1 0 0]); % median outer
                set(c(i+4*48), 'color', [1 0 0]); % box
                set(c(i+5*48), 'color', [1 0 0]); % whisker
            case 3 % target
                set(c(i+2*48), 'MarkerEdgeColor', [1 0.6 0.6]); % median inner
                set(c(i+3*48), 'MarkerEdgeColor', [1 0.6 0.6]); % median outer
                set(c(i+4*48), 'color', [1 0.6 0.6]); % box
                set(c(i+5*48), 'color', [1 0.6 0.6]); % whisker              
             
        end
    end
    
    % axis labels and properties
    xlim([-1 49]); ylim([0 350]);
    ylabel('Monthly rainfall (mm/month)')
    xlabels = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    for iMonth = 1:12, text(-2+4*iMonth, -25, xlabels{iMonth}, 'HorizontalAlignment', 'center'); end
    
    % % legend
    % text(27, 390, 'Legend', 'FontSize', 10)
    % text(28, 370, 'Blue: original (synthetic unperturbed)', 'FontSize', 8, 'color', [0 0 1])
    % text(28, 355, 'Light red: target', 'FontSize', 8, 'color', [1 0.6 0.6])
    % text(28, 340, 'Red: adopted', 'FontSize', 8, 'color', [1 0 0])
    % text(28, 325, 'Grey circles: outliers', 'FontSize', 8, 'color', [0.7 0.7 0.7])
    
    ThisFile = ['out\scratch\distributions_' RegionName '_forDeltaSeasonality' num2str(deltaSeasonality) '.png'];
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 17 30]); 
    % print('-painters','-depsc', ThisFile); 
    saveas(gcf, ThisFile, 'png');   
    close 
    
    dummy = -99.99;
    
end

function dummy = PlotSeasonalityPerturb2(orig, pert, RegionName, deltaSeasonality) % timeseries
        
    %% randomly choose a period to plot
    NumYears = 15; 
    start = 480; period = start:(start+12*NumYears-1); % 15 years
    
    close
    figure; hold on;
    plot(orig(period), 'b-');
    plot(pert(period), 'r-');
    
    % axis labelling and formatting
    ylabel('monthly rainfall (mm/month)'); ylim([0 200]);
    xlabel('year')
    XTickLabel_dbl = 1:NumYears; 
    for i = 1:NumYears, XTickLabel{i} = num2str(XTickLabel_dbl(i)); end
    set(gca, 'XTick', 1:12:(size(period, 2)+1), 'XTickLabel', XTickLabel)
    grid on
    
    legend('original (synthetic unperturbed)', 'adopted', 'Location', 'North')
    
    ThisFile = ['out\scratch\Timeseries_' RegionName '_forDeltaSeasonality' num2str(deltaSeasonality) '.png'];
    set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperPosition', [0 0 30 17]); 
    % print('-painters','-depsc', ThisFile); 
    saveas(gcf, ThisFile, 'png');   
    close 
    %%
    dummy = -99.99;
end