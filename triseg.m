function f = triseg(x,adjpars,data,init)
    
    %{ 

    This function contains the equations to calculate consistent initial
    conditions for the DAE. 

    Inputs: 
    x           - vector of states 
    pars        - vector of adjustable parameters 
    data        - input data structure with data and global parameters 
    init        - vector of initial conditions 

    Outputs: 
    f           - vector of equations for the root finder 

    %} 
    
    %% Unpack data structure 

    HR = data.HR; 
    fixpars = data.fixpars;

    %% Adjustable parameters 

    % Force scaling factors  
    k_pas_LV = adjpars(14);
    k_pas_RV = adjpars(15);
    k_act_LV = adjpars(16); 
    k_act_RV = adjpars(17);
    
    %% Fixed parameters
    
    % Sarcomere length parameters (um)
    Lsref   = fixpars(7);
    Lsc0    = fixpars(8); 
    Lse_iso = fixpars(9); 
    
    % Percentage of cardiac cycle 
    k_TS = fixpars(15); % Beginning of cardiac cycle to maximal systole  
    k_TR = fixpars(16); % Relaxation (maximal systole to baseline) 
    
    % Midwall reference surface area (cm^2)
    Amref_LV  = fixpars(19); 
    Amref_SEP = fixpars(20); 
    Amref_RV  = fixpars(21); 
    
    % Midwall volume (mL)
    Vw_LV  = fixpars(22); 
    Vw_SEP = fixpars(23); 
    Vw_RV  = fixpars(24); 

    %% Assign passive stress steepness parameter (dimensionless) 

    if isfield(data,'exvivo_on')
        gamma = adjpars(19); 
    else 
        gamma = fixpars(25); 
    end 
    
    if isfield(data,'eta_LVDD')
        eta_LVDD = data.eta_LVDD; 
    else 
        eta_LVDD = 1; 
    end 
    if isfield(data,'eta_RVDD')
        eta_RVDD = data.eta_RVDD; 
    else 
        eta_RVDD = 1; 
    end 

    %% Variables 
    
    % Undo log of states 
    x = exp(x); 

    % Displacements (cm) 
    xm_LV  = x(1); 
    xm_SEP = x(2); 
    xm_RV  = x(3); 
    ym     = x(4); 
    
    % Contractile element sarcomere length (um) 
    Lsc_LV  = x(5); 
    Lsc_SEP = x(6); 
    Lsc_RV  = x(7); 
    
    % Volumes (mL) 
    V_LV = init(8); 
    V_RV = init(9);
    
    %% Activation function 
    
    T = 60/HR; 
    TS = k_TS * T; 
    TR = k_TR * T; 
    
    t = 0; 
    tc_v = mod(t,T);
    if tc_v >= 0 && tc_v < TS 
        y_v = 0.5*(1 - cos(pi*tc_v/TS)); 
    elseif tc_v >= TS && tc_v < TR + TS 
        y_v = 0.5*(1 + cos(pi*(tc_v - TS)/TR)); 
    else
        y_v = 0; 
    end  
    
    %% Heart and sarcomere model 
    
    % Volume of spherical cap formed by midwall surface (m^3)
    Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
    Vm_SEP =  (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
    Vm_RV  =  (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 
    
    % Surface area of midwall surface (m^2) 
    Am_LV  = pi * (xm_LV^2  + ym^2);
    Am_SEP = pi * (xm_SEP^2 + ym^2); 
    Am_RV  = pi * (xm_RV^2  + ym^2); 
    
    % Curvature of midwall surface (m^(-1))
    Cm_LV  = -2 * xm_LV  / (xm_LV^2  + ym^2);
    Cm_SEP =  2 * xm_SEP / (xm_SEP^2 + ym^2);
    Cm_RV  =  2 * xm_RV  / (xm_RV^2  + ym^2);
    
    % Ratio of wall thickness to midwall radius of curvature (dimensionless)
    z_LV  = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV); 
    z_SEP = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP); 
    z_RV  = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV);
    
    % Myofiber strain (dimensionless)
    eps_LV  = 0.5 * log( Am_LV  / Amref_LV  ) - (1/12) * z_LV^2  - 0.019 * z_LV^4; 
    eps_SEP = 0.5 * log( Am_SEP / Amref_SEP ) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4; 
    eps_RV  = 0.5 * log( Am_RV  / Amref_RV  ) - (1/12) * z_RV^2  - 0.019 * z_RV^4; 
    
    % Sarcomere length (m)
    Ls_LV  = Lsref * exp(eps_LV); 
    Ls_SEP = Lsref * exp(eps_SEP); 
    Ls_RV  = Lsref * exp(eps_RV); 
    
    % Passive stress (kPa) 
    sigma_pas_LV  =  k_pas_LV * (Ls_LV - Lsc0)^(gamma * eta_LVDD); 
    sigma_pas_SEP =  k_pas_LV * (Ls_SEP - Lsc0)^(gamma*eta_LVDD); 
    sigma_pas_RV  =  k_pas_RV * (Ls_RV - Lsc0)^(gamma*eta_RVDD); 
    
    % Active stress (kPa)
    sigma_act_LV  = k_act_LV * y_v  * (Lsc_LV  - Lsc0) * (Ls_LV  - Lsc_LV)  / Lse_iso; 
    sigma_act_SEP = k_act_LV * y_v  * (Lsc_SEP - Lsc0) * (Ls_SEP - Lsc_SEP) / Lse_iso;
    sigma_act_RV  = k_act_RV * y_v  * (Lsc_RV  - Lsc0) * (Ls_RV  - Lsc_RV)  / Lse_iso;
    
    % Total stress (kPa)
    sigma_LV  = sigma_act_LV  + sigma_pas_LV; 
    sigma_SEP = sigma_act_SEP + sigma_pas_SEP; 
    sigma_RV  = sigma_act_RV  + sigma_pas_RV; 
    
    % Representative midwall tension (mmHg cm^(-2))
    Tm_LV  = (Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5); 
    Tm_SEP = (Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5); 
    Tm_RV  = (Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);
    
    % Axial midwall tension component 
    Tx_LV  = - Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2); 
    Tx_SEP =   Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2); 
    Tx_RV  =   Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2); 
    
    % Radial midwall tension component 
    Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2); 
    Ty_SEP = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2); 
    Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);
    
    %% System of equations 
    
    f1 = -V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV; 
    f2 = Tx_LV + Tx_SEP + Tx_RV;
    f3 = V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV;
    f4 = Ty_LV + Ty_SEP + Ty_RV; 
    
    f5 = (Ls_LV  - Lsc_LV)  /Lse_iso - 1;
    f6 = (Ls_SEP - Lsc_SEP) /Lse_iso - 1;
    f7 = (Ls_RV  - Lsc_RV)  /Lse_iso - 1;
    
    f = [f1; f2; f3; f4; f5; f6; f7]; 
end 

