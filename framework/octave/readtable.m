function outtable = readtable(filepath)
    
    % created 13/02/2021 by Keirnan Fowler, University of Melbourne
    
    % This is a *very basic* Octave readtable function to replicate readtable in Matlab.  
    % It is only intended to work with the data input files of the stochastic generation
    % framework example and may not be extendable to other contexts. 
    
    % suppress various pesky Octave warning messages (ensure you check
    % correct import though!)
    warning('off','all');
    
    % read the numerical data and column headers from the csv
    ThisData = importdata(filepath, ',', 1);
    
    % create separate vectors for each column
    NumCols = size(ThisData.colheaders, 2); namestring = [];
    for iCol = 1:NumCols
        ColName = ThisData.colheaders{iCol};
        eval([ColName ' = ThisData.data(:, iCol);']);
        namestring = [namestring ColName ', '];
    end
    
    % create table from variables
    eval(['outtable = table(' namestring(1:(end-2)) ');']);
    
    % turn warnings back on
    warning('default','all');
    
end