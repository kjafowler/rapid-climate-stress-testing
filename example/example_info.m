function info = example_info()
    
    % WAPABA parameter values for each region's representative catchment.
    % See Supplementary Material S6 for more information on 
    % how these parameter values were obtained.
    WapabaParSets.ParNames = {'a1', 'a2', 'B', 'Smax', 'K'};
    WapabaParSets.A  = [2.8468 1.7065 0.53284 174.1602   1.5096];
    WapabaParSets.B  = [1.9780 1.5908 0.52909 935.4290  24.0088];
    WapabaParSets.C  = [3.5389 1.8011 0.14338 231.0074   4.3214];
    WapabaParSets.D  = [2.8055 2.4473 0.38017 310.5108 207.6871];
    WapabaParSets.E  = [4.1223 1.6373 0.41727 199.4063  20.7406];
    WapabaParSets.F  = [3.3816 1.6207 0.32702 222.4516  16.3629];
    WapabaParSets.G  = [2.7690 1.5524 0.11801 173.1451  12.8342];
    
    % area of each subarea
    subarea = {'A', 'B', 'C', 'D', 'E', 'F', 'G'}';
    area_km2 = [3469 1291 3570 2423 3424 3240 26099]'; 
    
    % weighting to use when spatially averaging the low frequency component of precipitation
    % (ie weight by area but ignore region G - see paper, section 3.5)
    weighting_lowfreq = [3469 1291 3570 2423 3424 3240 0]';
    SubAreaDetails = table(subarea, area_km2, weighting_lowfreq); 
    
    % details of representative catchment for each subarea
    RepresentativeCatchment = {'405219', '405209', '405274', '405251', '406213', '407214', '407253'}';
    CA_km2 = [699.1 625.8 181.4 118.6 638.0 305.7 679.5]';
    identifier = {'A_405219', 'B_405209', 'C_405274', 'D_405251', 'E_406213', 'F_407214', 'G_407253'}';
    RepCatchDetails = table(subarea, RepresentativeCatchment, identifier, CA_km2); 
    
    % flow conversion factors - these are the multiplication factors 
    % applied to simulated flows from representative catchments, to obtain 
    % the reach inflows required by the river systems model.
    Eildoninflow = [4.10     0     0  3.50     0     0     0]';
    Reach1       = [   0  1.55  2.50     0     0     0     0]';
    Reach2       = [   0  0.80  7.00     0     0     0     0]';
    Reach3       = [   0     0  5.80     0  1.20     0  3.00]';
    Reach4       = [   0     0  7.00  2.50     0     0  2.00]';
    RepCatch = {'A_405219', 'B_405209', 'C_405274', 'D_405251', 'E_406213', 'F_407214', 'G_407253'}';
    FlowConversionFactors = table(RepCatch, Eildoninflow, Reach1, Reach2, Reach3, Reach4); 
    
    % other parameters
    pars.NumSubareas = size(RepCatchDetails, 1); 
    pars.WaterYearStart_clim = 'January'; % basis of water years for climate ('January' means calendar years are used)
    pars.WaterYearStart_flow = 'March';   % basis of water years for streamflow
    pars.StochRepLen_yrs  =  3000;        % length of stochastic replicates, in years
    pars.Streamflow_BoxCox_Lambda = 0.79; % used during perturbation of the rainfall-runoff relationship (Section 2.4.5 / Table 1 / Section 3.7 / Supp Mat Section S7)
    
    % store all the above in structure 'info'
    info = struct('WapabaParSets', WapabaParSets, 'SubAreaDetails', SubAreaDetails, 'RepCatchDetails', RepCatchDetails, 'FlowConversionFactors', FlowConversionFactors, 'pars', pars);
    info.SubareaList = subarea;
    
end
