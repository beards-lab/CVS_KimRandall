function y = activation(t,data)

%% Unpack the data structure 
    
    HR      = data.HR; 
    fixpars = data.fixpars; 

%% Fixed parameters 

    % Percentage of cardiac cycle 
    k_TS = fixpars(15); % Beginning of cardiac cycle to maximal systole  
    k_TR = fixpars(16); % Relaxation (maximal systole to baseline 

%% Activation function
    
    % Heart period (s) 
    T = 60/HR; 
    
    % Time to maximal systole 
    TS_v = k_TS * T; 
    
    % Time from maximal systole to relaxation 
    TR_v = k_TR * T; 
    
    % Ventricular activation function 
    tc_v = mod(t,T);
    if tc_v >= 0 && tc_v < TS_v 
        y = 0.5*(1 - cos(pi*tc_v/TS_v)); 
    elseif tc_v >= TS_v && tc_v < TR_v + TS_v 
        y = 0.5*(1 + cos(pi*(tc_v - TS_v)/TR_v)); 
    else
        y = 0; 
    end 

end 