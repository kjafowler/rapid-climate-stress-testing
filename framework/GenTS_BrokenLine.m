function TS_BrokenLine = GenTS_BrokenLine(param_timescale, param_amplitude, TS_rand, pars)

    % generate the stochastic timeseries for the low frequency component, 
    % using the Broken Line process with the parameter values specified
    
    n = pars.StochRepLen_yrs; 
    
    % generate the start and end of the meta-timesteps...
    i = 1; % segment counter
    j = 1; % year counter
    ThisVal = norminv(TS_rand(i+1), 0, param_amplitude); 
    est_size = floor(n/param_timescale); 
    st(1:est_size)=NaN; en(1:est_size)=NaN; lookup_ij(1:n) = NaN; % saves time
    
    TS_end = 1; 
    while TS_end < n
        i = i+1;
        TS_start = TS_end;
        TS_end = TS_start + param_timescale;
        st(i) = TS_start; 
        en(i) = TS_end; 
        
        % ...along with the values at the transition points
        val_st(i) = ThisVal; % the end value of the previous metastep is the start value of this one
        ThisVal = norminv(TS_rand(i+1), 0, param_amplitude); % use TS_rand: pre-generated random numbers
        val_en(i) = ThisVal; 
        
        % and create a table so that we can tell which timesteps are in
        % this meta-timestep
        while j<=TS_end
            lookup_ij(j) = i;
            j = j+1; 
        end
    end
    
    % generate the annual series
    TS_BrokenLine(1:n) = NaN; 
    for j = 1:n
        i = lookup_ij(j); % max(find(j>=st));
        x1 = st(i);
        x2 = en(i);
        y1 = val_st(i);
        y2 = val_en(i);
        x  = j;
        assert(abs((x2-x1) - param_timescale)<0.000001);
        TS_BrokenLine(j)= y1+ ((x-x1)/(x2-x1)) * (y2-y1);
    end    

end