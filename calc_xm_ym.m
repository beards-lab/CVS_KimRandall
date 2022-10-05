function f = calc_xm_ym(x,Lsref,Vw,Amref,SL,V_LV,flag)

    %{ 

    This function is used in the calculation of the initial estimates for
    the xm and ym values and the midwall surface area reference parameter
    for the RV (Amref_RV). 

    Inputs: 
    x       - vector of states 
    Lsref   - reference sarcomere length parameter 
    Vw      - vector of wall volume parameters 
    Amref   - vector of midwall surface area reference parameters 
    SL      - sarcomere length 
    V_LV    - left ventricular volume (commonly the end-diastolic volume) 
    flag    - flag to determine whether calculation of end-diastole or
    end-systole 

    %}

    %% Parameters 
    
    % Midwall volumes (mL)
    Vw_LV  = Vw(1); 
    Vw_SEP = Vw(2); 
    Vw_RV  = Vw(3); 
    
    % Midwall reference surface areas (cm^2)
    Amref_LV  = Amref(1); 
    Amref_SEP = Amref(2); 
    
    % If in end-systole, use assigned Amref_RV value 
    if flag == 1 % 0 = end-diastole ; 1 = end-systole 
        Amref_RV = Amref(3); 
    end 
    
    %% States
    
    x = exp(x); 
    
    % Displacements (cm) 
    xm_LV    = x(1); 
    xm_SEP   = x(2); 
    xm_RV    = x(3);
    ym       = x(4); 
    
    % If in end-diastole, treat Amref_RV as a state 
    if flag == 0 
        Amref_RV = x(5); 
    end 
    
    %% Equations
    
    % Midwall surface volume (mL)
    Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
    Vm_SEP = (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
    Vm_RV  = (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 
    
    % Midwall surface area (cm^2)
    Am_LV  = pi * (xm_LV^2  + ym^2);
    Am_SEP = pi * (xm_SEP^2  + ym^2);
    Am_RV  = pi * (xm_RV^2  + ym^2);
    
    % Midwall curvature (cm^(-1))
    Cm_LV  = -2 * xm_LV  / (xm_LV^2  + ym^2);
    Cm_SEP = -2 * xm_SEP / (xm_SEP^2  + ym^2);
    Cm_RV  = -2 * xm_RV  / (xm_RV^2  + ym^2);
    
    % Midwall ratio (dimensionless) 
    z_LV   = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV); 
    z_SEP  = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP); 
    z_RV   = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV); 
    
    % Strain (dimensionless) 
    eps_LV  = 0.5 * log( Am_LV  / Amref_LV  ) - (1/12) * z_LV^2  - 0.019 * z_LV^4; 
    eps_SEP = 0.5 * log( Am_SEP / Amref_SEP ) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4; 
    eps_RV  = 0.5 * log( Am_RV  / Amref_RV  ) - (1/12) * z_RV^2  - 0.019 * z_RV^4; 
    
    % Instantaneous sarcomere length (um) 
    Ls_LV  = Lsref * exp(eps_LV); 
    Ls_SEP = Lsref * exp(eps_SEP); 
    Ls_RV  = Lsref * exp(eps_RV); 
    
    V_RV = V_LV; 
    
    %% Outputs
    
    % Balance volumes 
    f(1) = abs(-V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV); 
    f(2) = abs(V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV);
    
    % Balance sarcomere lengths in end-diastole 
    if flag == 0 
        f(3) = abs(Ls_LV  - SL);
        f(4) = abs(Ls_SEP - SL); 
        f(5) = abs(Ls_RV  - SL); 
    end 
    
    end 