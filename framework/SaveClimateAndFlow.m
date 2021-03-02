
function SaveClimateAndFlow(i, DataOut, deltaP, deltaT, deltaLowFreqP, deltaSeasonality, deltaRRrship_space, info)          
    
    if ~info.isOctave
        % get the data relevant to this deltaRRrship
        deltaRRrship = deltaRRrship_space(i);
        ThisData = DataOut.(['deltaRRrship' num2str(i)]); 

        % separate out the items to be saved from the input structure
        TS_Q_forModel = ThisData.TS_Q_ReachInflows;
        TS_P          = ThisData.TS_P; 
        TS_T          = ThisData.TS_T; 
        TS_PET        = ThisData.TS_PET; 

        % remember the perturbations
        perturbations = struct('deltaP', deltaP, 'deltaT', deltaT, 'deltaLowFreqP', deltaLowFreqP, 'deltaSeasonality', deltaSeasonality, 'deltaRRrship', deltaRRrship); 

        % create name of file
        OutFile_path = ['out\out_data\' 'StochasticClimate_' ...
            num2str(deltaP, '%+3.2f') '_' ...
            num2str(deltaT, '%+3.2f') '_' ...
            num2str(deltaLowFreqP, '%+3.3f') '_' ...
            num2str(deltaSeasonality, '%+3.2f') '_' ...
            num2str(deltaRRrship, '%+3.2f') '.mat'];

        % now save
        save(OutFile_path, 'TS_Q_forModel', 'TS_P', 'TS_T', 'TS_PET', 'perturbations');
    else
        warning('Note: Octave is not currently able to save tables to .mat files.  Users of Octave should code their own procedure to archive the output.');
    end
end

