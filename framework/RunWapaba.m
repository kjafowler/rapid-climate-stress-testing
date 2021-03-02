function [SimFlows, SimVars] = RunWapaba(ParVals, rain, PET, month)

    % Created 10/05/2018 by Connor McCutcheon and Keirnan Fowler
    
    % Purpose: run Wapaba model with specified timeseries inputs and model parameters
    
    % Model Parameters
    a1 = ParVals.a1;
    a2 = ParVals.a2;
    B = ParVals.B;
    Smax = ParVals.Smax;
    K = ParVals.K;
   
    % initialise size of arrays
    NumTimesteps = size(rain, 1); 
    SimFlows         (1:NumTimesteps) = -99.99; 
    T                (1:NumTimesteps) = -99.99; 
    SimVars.Flows.Qb (1:NumTimesteps) = -99.99; 
    SimVars.Flows.Qs (1:NumTimesteps) = -99.99; 
    SimVars.States.S (1:NumTimesteps) = -99.99; 
    SimVars.States.G (1:NumTimesteps) = -99.99; 
    SimVars.Flux.ET  (1:NumTimesteps) = -99.99; 
    SimVars.Flux.R   (1:NumTimesteps) = -99.99; 
    
    % Step 1: initialise model states
    S = 0.01*Smax; 
    G = 0;         
    
    % Step 2: loop to run model for each timestep
    for iTS = 1:NumTimesteps
        
       P = rain(iTS);

      % Define T value for each month
       if month (iTS) == 1 || 3 || 5 || 7 || 8 || 10 || 12
           T (iTS) = 31;
            elseif month (iTS) == 2
                T (iTS) = 28;
            else 
                T (iTS) = 30;
       end    
  
       %Catchment water consumption potential 
        X0 = PET(iTS) + (Smax - S);
       
       %Consumption curve #1: Supply = rainfall, Demand = Catchment water
       %consumption potential
        F1 = 1+P/X0-(1+(P/X0).^a1).^(1/a1);
       
       %Catchment water consumption
        X = X0 * F1;
       
       %Catchment water yield
        Y = max(P - X,0);
       
       %Total water available for evapotranspiration
        W = S + X;
       
       %Consumption curve #2: Supply = water available for ET, Demand = PET
        F2 = 1+W/PET(iTS)-(1+(W/PET(iTS)).^a2).^(1/a2);
       
       %Actual evapotranspiration
        ET = F2 * PET(iTS);
       
       %Recharge to groundwater store
        R = max(B * Y,0);
       
       %Surface Runoff
        Qs = max(Y - R,0);
       
       %Baseflow calculation
        Z = 1-exp(-T(iTS)/K);
        Qb1 = min(G*Z+R*(1-(K/T(iTS))*Z),G);
        Qb= max(Qb1,0);
       
       %Groundwater store
        G = max(G+R-Qb,0);
       
       %Soil moisture storage
        S1 = max(W-ET,0);
        S= min(S1,Smax);
       
        %Accounting for soil moisture storage exceedance
       if S1>Smax
           Qs=Qs+(S1-Smax);
       end 
        
        %Total Runoff
         Q = Qs + Qb;
        
        % Assign SimFlows and SimVars
        SimFlows(iTS) = Q; 
        SimVars.States.S(iTS) = S; 
        
        % some optional tasks
        RememberAdditionalVariables = false; % option to turn on if user requires
        if RememberAdditionalVariables
            SimVars.Flows.Qb (iTS) = Qb;
            SimVars.Flows.Qs (iTS) = Qs;        
            SimVars.States.G(iTS) = G; 
            SimVars.Flux.ET(iTS) = ET;
            SimVars.Flux.R(iTS) = R;
        end
        
        % optional water balance checking
        WaterBalanceChecking = false; % option to turn on if user requires
        if WaterBalanceChecking
            A=iTS-1;
            if A==0
                SystemBal=0;
            else
              % SystemBal = round(P + SimVars.States.S(iTS-1) - SimVars.States.S(iTS) + SimVars.States.G(iTS-1) - SimVars.States.G(iTS) - SimVars.Flux.ET(iTS) - SimFlows(iTS),5);
                SystemBal = round(P + LastTS.S                - S                     + LastTS.G                - G                     - ET                   - Q,5);
            end
            SimVars.Balance(iTS) = SystemBal;
        end
        
        if isreal(S) && isreal(G) && isreal(Q)
        else
            if ~exist('ComplexNumberDetected')
                ParVals
                error(['Complex values in simulated Q for timestep ' num2str(iTS) ' for above parameter set.']);
                ComplexNumberDetected = true;
            end
        end
    end
  
end