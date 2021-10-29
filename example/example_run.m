
% Created January 2021 by Keirnan Fowler, University of Melbourne, fowler.k@unimelb.edu.au

% Integrated framework for rapid climate stress testing on a monthly timestep
% by Keirnan Fowler, Natasha Ballis, Avril Horne, Andrew John, Rory Nathan and Murray Peel 

% Licence: CC BY 3.0 - see https://creativecommons.org/licenses/by/3.0/au/

% Use this file to run the example given in the paper.  Specifically, this 
% file can be used to demonstrate for the following two cases: 
% - case 1: only generate the stochastic data required to create the 
%           plots in the paper, and no more.  Then do the plotting.
% - case 2: generate a full dataset of stochastic data covering the entire
%           stress testing space, as defined by the axis limits and 
%           gradations specified in Table # in the paper

% This file uses the code contained in folder '\framework', producing 
% stochastic output saved in folder '\example\out' (in .mat format), then 
% using code from example_plotting.m to create plots.

clear; close all; 

% add path to directories with framework code
addpath('..\framework\');
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; % test if we are in Octave
if isOctave, addpath(genpath('..\framework\octave\')); pkg load statistics; end

% start timing the code
tic;

%% load inputs including historic data

% climatic data
precip_monthly = readtable('hist\precip_monthly.csv'); 
T_monthly = readtable('hist\Tmax_monthly.csv'); 
PET_monthly = readtable('hist\PET_monthly.csv'); 

% flow data - note, this is not used directly by the framework but is used
% to produce plots that appear in the paper.  The unit of flow varies with
% context.
flow_monthly.RepresentativeCatchments_mm = readtable('hist\flow_monthly_representative_catchments.csv'); 
flow_monthly.ReachInflows_ML = readtable('hist\flow_monthly_reach_inflows.csv'); 

% store all this data in a single structure
HistoricData = struct('precip_monthly', precip_monthly, 'T_monthly', T_monthly, 'PET_monthly', PET_monthly, 'flow_monthly', flow_monthly);

info = example_info();  % specify information and settings, including the following: 
                        % info.WapabaParSets          Parameter sets for pre-calibrated WAPABA rainfall runoff model - see Supplementary Material ##
                        % info.SubareaDetails         Subarea names, areas, and weightings to use in aggregation of Intrinsic Mode Functions (IMFs) (see function AggregateIMFs)
                        % info.RepCatchDetails        Each subarea has a "representative" catchment.  This table stores the characteristics of these catchments.  See paper Section ## and Supplementary Material ##.  
                        % info.FlowConversionFactors  Each reach of the Goulburn River receives inflows from multiple subareas.  This table of factors describes how.  See paper Section ## and Supplementary Material ##.   
                        % info.pars                   Other general settings relevant to stochastic climate generation and the way the data is used
                    
% aggregate all data to annual
% note, the results are appended to the existing data structure 'HistoricData'
HistoricData = AggregateToAnnual(HistoricData, info.pars); 

info.isOctave = isOctave; 
clearvars -except HistoricData info

%% Split precipitation into high and low frequency using Empirical Mode Decomposition
% use the CEEMDAN algorithm (Complete Ensemble Empirical Mode Decomposition with Adaptive Noise) 
% note, the results are appended to the existing data structure 'HistoricData'
disp('Finished loading data.  Starting empirical mode decomposition using CEEMDAN...')

% set the threshold that determines which IMFs are 'low' and 'high' (see paper, Section ##) 
info.LowHighThresh = 2; % IMFs 1 and 2 are high frequency, everything else is low

% run CEEMDAN
HistoricData = split_high_low_using_CEEMDAN(HistoricData, info);

%% Conduct pre-analysis to inform stochastic generation of the low-frequency component
disp('Done.  Starting pre-analysis for low-frequency component.');
info.LowFreq_PreAnalysis_Outputs = LowFreq_PreAnalysis(HistoricData, info); 

timing.init = toc; tic;

%% case 1: just the plots in the paper
%  only generate the stochastic data required to create the 
%  plots in the paper, and no more.  Then do the plotting.
disp('Done.  Starting stochastic generation...');

for i = 1:12
    
    % only test eleven different combinations of stressors (ie. the ones in Figure 10)
    switch i
        case  1, name = 'BaseCase';           deltaP = 0.00; deltaT =  0; deltaLowFreqP =      0; deltaSeasonality =     0; deltaRRrship =     0; 
        case  2, name = 'deltaP_down';        deltaP = -0.3; deltaT =  0; deltaLowFreqP =      0; deltaSeasonality =     0; deltaRRrship =     0; 
        case  3, name = 'deltaP_up';          deltaP = +0.1; deltaT =  0; deltaLowFreqP =      0; deltaSeasonality =     0; deltaRRrship =     0; 
        case  4, name = 'deltaT_up';          deltaP =  0.0; deltaT = +2; deltaLowFreqP =      0; deltaSeasonality =     0; deltaRRrship =     0; 
        case  5, name = 'deltaLowFreqP_up';   deltaP =  0.0; deltaT =  0; deltaLowFreqP = + 0.06; deltaSeasonality =     0; deltaRRrship =     0;         
        case  6, name = 'deltaLowFreqP_down'; deltaP =  0.0; deltaT =  0; deltaLowFreqP = -0.015; deltaSeasonality =     0; deltaRRrship =     0;         
        case  7, name = 'deltaSeas_up';       deltaP =  0.0; deltaT =  0; deltaLowFreqP =      0; deltaSeasonality = +0.12; deltaRRrship =     0;         
        case  8, name = 'deltaSeas_down';     deltaP =  0.0; deltaT =  0; deltaLowFreqP =      0; deltaSeasonality = -0.06; deltaRRrship =     0;         
        case  9, name = 'deltaRRrel_up';      deltaP =  0.0; deltaT =  0; deltaLowFreqP =      0; deltaSeasonality =     0; deltaRRrship =   -25;         
        case 10, name = 'deltaRRrel_down';    deltaP =  0.0; deltaT =  0; deltaLowFreqP =      0; deltaSeasonality =     0; deltaRRrship = +6.25;         
        case 11, name = 'combination';        deltaP = -0.1; deltaT =  1; deltaLowFreqP = +0.015; deltaSeasonality = +0.03; deltaRRrship = -6.25; 
        case 12, name = 'combination2';       deltaP = -0.2; deltaT =  2; deltaLowFreqP = +0.030; deltaSeasonality = +0.06; deltaRRrship = -12.5; 
    end
    
    % generate stochastic perturbed data for this combination of stressors
    DataOut = GetStochPertData(deltaP, deltaT, deltaLowFreqP, deltaSeasonality, deltaRRrship, HistoricData, info);
    
    % store in a common data structure
    StochPertData.(name) = DataOut; 
    disp(['Case 1 stochastic generation: ' num2str(i) ' of 11 done'])
    
end
timing.case1_stochgen = toc; tic;

% create plots
example_plotting(StochPertData, HistoricData, info); 
timing.plotting = toc; tic;

%% case 2: full stochastic dataset
%  generate a full dataset of stochastic data covering the entire
%  stress testing space, as defined by the axis limits and 
%  gradations specified in Table # in the paper.  Whereas the code above
%  placed all the outputs into a single structure (StochPertData), here
%  there are too many so we save them as individual .mat files.  

% specify axis gradations
deltaP_space           = [-0.40 -0.35 -0.30 -0.25 -0.20 -0.15 -0.10 -0.05 0 .05 .10 .15]; % rainfall proportional change
deltaT_space           = [0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0];                             % additional degrees of warming
deltaLowFreqP_space    = [-0.03 -0.015 0 +0.015 +0.03 +0.045 +0.06 +0.075];               % changes n Hurst Coefficient
deltaSeasonality_space = [-0.06 -0.03 0.00 +0.03 +0.06 +0.09 +0.12 +0.15];                % changes in seasonality
deltaRRrship_space     = [-50 -43.75 -37.5 -31.25 -25 -18.75 -12.5 -6.25 0 6.25 12.5];    % shift in rainfall runoff relationship

% generate data
for deltaP = deltaP_space 
    for deltaT = deltaT_space
        for deltaLowFreqP = deltaLowFreqP_space
            for deltaSeasonality = deltaSeasonality_space
                
                DesiredNumberOfTestRuns = 1000;
                if rand() < (DesiredNumberOfTestRuns / (12*9*7*8*11))
                    
                    % run the stochastic data routines
                    DataOut = GetStochPertData(deltaP, deltaT, deltaLowFreqP, deltaSeasonality, deltaRRrship_space, HistoricData, info);

                    % NOTE: 
                    % The following notes relate to function GetStochPertData, which is run both here and above at line 82.  
                    % At line 82, the input 'deltaRRrship' was a single value, whereas here it is a vector.  In either case
                    % the stochastic climate generation is run only once.  This reflects that the stochastic climate generation
                    % tasks do not depend on deltaRRrship and thus it is a waste of effort to rerun them for each new value.  
                    % In contrast, the flow perturbation code is run either once (in the case of line 82) or multiple times 
                    % in the case of a vector (once for each value in the vector).   

                    % save to file
                    for i = 1:size(deltaRRrship_space, 2)
                        SaveClimateAndFlow(i, DataOut, deltaP, deltaT, deltaLowFreqP, deltaSeasonality, deltaRRrship_space, info); 
                    end
                    % *note, to save disc space the user could alter this to only save one copy of the climate inputs at this point, and then
                    %  a separate file for each deltaRRrship.  
                end
            end
        end
    end
end
timing.case2_stochgen = toc; 

%% report back on timing
disp('All done.  Run times were:')
disp('Initialisation including data loading and low frequency pre-analysis: ')
disp([num2str(timing.init) ' seconds. ']); 
disp('Stochastic generation for case 1')
disp([num2str(timing.case1_stochgen) ' seconds. ']); 
disp('Plotting for case 1')
disp([num2str(timing.plotting) ' seconds. ']); 
disp('Stochastic generation for case 2')
disp([num2str(timing.case2_stochgen) ' seconds. ']); 
