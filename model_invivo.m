function [dxdt, outputs] = model_invivo(t,x,adjpars,data) 

    %{ 

    This function contains the right-hand side of the in vivo model. 

    Inputs: 
    t       - time 
    x       - states 
    pars    - vector of adjustable parameter values 
    data    - input data structure with data and global parameters 

    Outputs: 
    dxdt    - vector of solved ODEs
    outputs - vector of solved auxiliary equations 

    %}

    %% Unpack data structure 

    HR = data.HR; 
    fixpars = data.fixpars;
    
    %% Adjustable parameters
    
    % Compliance (mL mmHg^(-1))
    C_SA = adjpars(1); 
    C_SV = adjpars(2); 
    C_PA = adjpars(3); 
    C_PV = adjpars(4); 
    
    % Resistance (mmHg s mL^(-1))
    R_SA  = adjpars(5); 
    R_tSA = adjpars(6); 
    R_PA  = adjpars(7); 
    
    R_m_vlv = adjpars(8); 
    R_a_vlv = adjpars(9); 
    R_t_vlv = adjpars(10); 
    R_p_vlv = adjpars(11);  
    
    % Pericardium parameters 
    Vh0 = adjpars(12); % mL 
    s   = adjpars(13);
    
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
    
    % Sarcomere length shortening velocity (um s^(-1))
    v_max = fixpars(10); 
    
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
    
    % Passive stress steepness parameter (dimensionless) 
    gamma = fixpars(25); 

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
    
    % Axial distance of midwall junction (cm)
    xm_LV  = x(1); 
    xm_SEP = x(2); 
    xm_RV  = x(3);
    
    % Radial distance of midwall junction (cm)
    ym = x(4); 
    
    % Contractile element length (um)
    Lsc_LV  = x(5); 
    Lsc_SEP = x(6); 
    Lsc_RV  = x(7); 
    
    % Volumes (mL) 
    V_LV = x(8); 
    V_RV = x(9);
    V_SA = x(10); 
    V_SV = x(11);
    V_PA = x(12); 
    V_PV = x(13); 
    
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
        y_v = 0.5*(1 - cos(pi*tc_v/TS_v)); 
    elseif tc_v >= TS_v && tc_v < TR_v + TS_v 
        y_v = 0.5*(1 + cos(pi*(tc_v - TS_v)/TR_v)); 
    else
        y_v = 0; 
    end 
    
    %% Heart model
    
    % Volume of spherical cap formed by midwall surface (mL)
    Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
    Vm_SEP = (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
    Vm_RV  = (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 
    
    % Surface area of midwall surface (cm^2) 
    Am_LV  = pi * (xm_LV^2  + ym^2);
    Am_SEP = pi * (xm_SEP^2 + ym^2); 
    Am_RV  = pi * (xm_RV^2  + ym^2); 
    
    % Curvature of midwall surface (cm^(-1))
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
    
    % Sarcomere length (um)
    Ls_LV  = Lsref * exp(eps_LV); 
    Ls_SEP = Lsref * exp(eps_SEP); 
    Ls_RV  = Lsref * exp(eps_RV); 
    
    % Passive stress  
    sigma_pas_LV  =  k_pas_LV * (Ls_LV - Lsc0)^(gamma*eta_LVDD); 
    sigma_pas_SEP =  k_pas_LV * (Ls_SEP - Lsc0)^(gamma*eta_LVDD); 
    sigma_pas_RV  =  k_pas_RV * (Ls_RV - Lsc0)^(gamma*eta_RVDD); 
    
    % Active stress 
    sigma_act_LV  = k_act_LV * y_v  * (Lsc_LV  - Lsc0) * (Ls_LV  - Lsc_LV)  / Lse_iso; 
    sigma_act_SEP = k_act_LV * y_v  * (Lsc_SEP - Lsc0) * (Ls_SEP - Lsc_SEP) / Lse_iso;
    sigma_act_RV  = k_act_RV * y_v  * (Lsc_RV  - Lsc0) * (Ls_RV  - Lsc_RV)  / Lse_iso;
    
    % Total stress 
    sigma_LV  = sigma_act_LV  + sigma_pas_LV; 
    sigma_SEP = sigma_act_SEP + sigma_pas_SEP; 
    sigma_RV  = sigma_act_RV  + sigma_pas_RV; 
    
    % Representative midwall tension 
    Tm_LV  = (Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5); 
    Tm_SEP = (Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5); 
    Tm_RV  = (Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);
    
    % Axial midwall tension component 
    Tx_LV  = - Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2); 
    Tx_SEP = Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2); 
    Tx_RV  = Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2); 
    
    % Radial midwall tension component 
    Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2); 
    Ty_SEP = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2); 
    Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);
    
    % Ventricular pressure (kPa)
    ptrans_LV = 2 * Tx_LV / ym; 
    ptrans_RV = 2 * Tx_RV / ym; 
    P_LV      = -ptrans_LV; 
    P_RV      = ptrans_RV; 
    
    %% Lumped circulatory model 

    % Pericardial pressure (mmHg) 
    Vh = V_LV + V_RV; 
    P_peri = exp(s*(Vh/Vh0 - 1)); % convert mmHg to kPa 
    P_LV = P_peri + P_LV; 
    P_RV = P_peri + P_RV; 
    
    % Pressure (mmHg)
    P_SV = V_SV / C_SV; 
    P_PA = V_PA / C_PA; 
    P_PV = V_PV / C_PV;  
    
    % When aortic valve is closed 
    Q_a_vlv = 0; 
    P_SA = (R_SA*V_SA + C_SA*P_SV*R_tSA + C_SA*Q_a_vlv*R_SA*R_tSA)/(C_SA*(R_SA + R_tSA));
    Q_SA = (V_SA - C_SA*P_SV + C_SA*Q_a_vlv*R_tSA)/(C_SA*(R_SA + R_tSA));
    
    % When aortic valve is open 
    if (P_SA < P_LV) * (V_LV > 0)    
        Q_a_vlv = -(R_SA*V_SA - C_SA*P_LV*R_SA - C_SA*P_LV*R_tSA + C_SA*P_SV*R_tSA)/(C_SA*(R_SA*R_a_vlv + R_SA*R_tSA + R_a_vlv*R_tSA));
        Q_SA = (R_a_vlv*V_SA - C_SA*P_SV*R_a_vlv + C_SA*P_LV*R_tSA - C_SA*P_SV*R_tSA)/(C_SA*(R_SA*R_a_vlv + R_SA*R_tSA + R_a_vlv*R_tSA));
        P_SA = (R_SA*R_a_vlv*V_SA + C_SA*P_LV*R_SA*R_tSA + C_SA*P_SV*R_a_vlv*R_tSA)/(C_SA*(R_SA*R_a_vlv + R_SA*R_tSA + R_a_vlv*R_tSA));
    end
    

    % Flows (mL s^(-1))
    Q_m_vlv = max((P_PV - P_LV) / R_m_vlv, 0); 
    Q_t_vlv = max((P_SV - P_RV) / R_t_vlv, 0); 
    Q_p_vlv = max((P_RV - P_PA) / R_p_vlv, 0); 
    Q_PA = (P_PA - P_PV) / R_PA; 

    %% ODEs
    
    % 1 - 4
    dxm_LV  = -V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV; 
    dxm_SEP = Tx_LV + Tx_SEP + Tx_RV;
    dxm_RV  = V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV;
    dym     = Ty_LV + Ty_SEP + Ty_RV; 
    
    % 5 - 7
    dLsc_LV  = ((Ls_LV  - Lsc_LV)  / Lse_iso - 1) * v_max;
    dLsc_SEP = ((Ls_SEP - Lsc_SEP) / Lse_iso - 1) * v_max;
    dLsc_RV  = ((Ls_RV  - Lsc_RV)  / Lse_iso - 1) * v_max;
    
    % 8 - 14
    dV_LV = Q_m_vlv - Q_a_vlv; 
    dV_SA = Q_a_vlv - Q_SA; 
    dV_SV = Q_SA    - Q_t_vlv; 
    dV_RV = Q_t_vlv - Q_p_vlv; 
    dV_PA = Q_p_vlv - Q_PA; 
    dV_PV = Q_PA    - Q_m_vlv; 
    
    dxdt = [dxm_LV; dxm_SEP; dxm_RV; dym;
        dLsc_LV; dLsc_SEP; dLsc_RV; 
        dV_LV; dV_RV; dV_SA; dV_SV; dV_PA; dV_PV; 
        ]; 
    
    outputs = [P_LV; P_SA; P_SV; P_RV; P_PA; P_PV;          % 1-6
        Vm_LV; Vm_SEP; Vm_RV;                               % 7-9
        Am_LV; Am_SEP; Am_RV;                               % 10-12
        Cm_LV; Cm_SEP; Cm_RV;                               % 13-14
        eps_LV; eps_SEP; eps_RV;                            % 15-18
        sigma_pas_LV; sigma_pas_SEP; sigma_pas_RV;          % 19-21
        sigma_act_LV; sigma_act_SEP; sigma_act_RV;          % 22-24
        sigma_LV; sigma_SEP; sigma_RV;                      % 25-27
        Q_m_vlv; Q_a_vlv; Q_t_vlv; Q_p_vlv;                 % 28-31
        Q_SA; Q_PA;                                         % 32-33
        Tm_LV; Tm_SEP; Tm_RV;                               % 34-36
        y_v;                                                % 37
        ];

end 