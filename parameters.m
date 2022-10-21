function [adjpars,UB,LB,data] = parameters(data)

    %{ 
    Assignment and/or nominal calculation of all model parameters. 
    
    Inputs: 
    data        - input data structure with data and global parameters 

    Outputs: 
    adjpars     - vector of adjustable parameters 
    UB          - vector of parameter upper bounds 
    LB          - vector of parameter lower bounds 
    data        - output data structure with new field assignments 

    %} 

    %% Unpack data structure
    
    % Blood pressures (mmHg)
    SPbar = data.SPbar; 
    DPbar = data.DPbar; 
    
    % Total blood volume (mL) 
    Vtot  = data.Vtot; 
    
    % Cardiac output (mL s^(-1))
    CO    = data.CO;
    
    % End-diastolic and end-systolic pressures (mmHg) and volumes (mL) 
    ESP_LV = data.ESP_LV;   ESP_RV = data.ESP_RV; 
    ESV_LV = data.ESV_LV;    
    EDP_LV = data.EDP_LV;   EDP_RV = data.EDP_RV;
    EDV_LV = data.EDV_LV;   EDV_RV = data.EDV_RV;
    
    %% Volumes
    
    % Blood volume distribution values - sum total = 1.0 
    d_LV = .025;             d_RV = .025;
    d_SA = .15;              d_PA = .05; 
    d_SV = .55;              d_PV = .2;

    % Total chamber volumes 
    V_LV_0 = d_LV*Vtot;      V_RV_0 = d_RV*Vtot; 
    V_SA_0 = d_SA*Vtot;      V_PA_0 = d_PA*Vtot; 
    V_SV_0 = d_SV*Vtot;      V_PV_0 = d_PV*Vtot;
    
    % Unstressed volume percentages for arterial and venous compartments 
    bvd_SA = data.bvd_SA;   bvd_PA = data.bvd_PA; 
    bvd_SV = data.bvd_SV;   bvd_PV = data.bvd_PV; 

    % Unstressed chamber volumes
    V_SA_u = V_SA_0*bvd_SA;   V_PA_u = V_PA_0*bvd_PA; 
    V_SV_u = V_SV_0*bvd_SV;   V_PV_u = V_PV_0*bvd_PV; 
    
    % Stressed chamber volumes
    V_SA_s = V_SA_0 - V_SA_u;  V_PA_s = V_PA_0 - V_PA_u;  
    V_SV_s = V_SV_0 - V_SV_u;  V_PV_s = V_PV_0 - V_PV_u; 
    
    %% Pressures
    
    % General max and min chamber pressures from Boron
    % LV            % SA                % SV
    P_LV_M = 125;   P_SA_M   = 120;     P_SV_M = 6;
                    P_SA_bar = 95;      P_SV_bar = 4; 
    P_LV_m = 1;     P_SA_m   = 80;      P_SV_m = 1;
                                                
    % RV            % PA                % PV
    P_RV_M = 26;    P_PA_M = 20;        P_PV_M = 8; 
                    P_PA_bar = 15;      P_PV_bar = 4; 
    P_RV_m = 1;     P_PA_m = 10;        P_PV_m = 2; 
    
    % Scale pressures to max and min SA pressures 
    q_LV_M = P_LV_M/P_SA_M;    q_LV_m = P_LV_m/P_SA_m; 
    q_RV_M = P_RV_M/P_SA_M;    q_RV_m = P_RV_m/P_SA_m; 
    q_SA_M = P_SA_M/P_SA_M;    q_SA_m = P_SA_m/P_SA_m; 
    q_SV_M = P_SV_M/P_SA_M;    q_SV_m = P_SV_m/P_SA_m; 
    q_PA_M = P_PA_M/P_SA_M;    q_PA_m = P_PA_m/P_SA_m; 
    q_PV_M = P_PV_M/P_SA_M;    q_PV_m = P_PV_m/P_SA_m; 
    
    % Scale chamber pressures to the mean SBP and DBP from the data     
    P_LV_M = q_LV_M*SPbar;    P_LV_m = q_LV_m*DPbar;       
    P_RV_M = q_RV_M*SPbar;    P_RV_m = q_RV_m*DPbar;   
    P_SA_M = q_SA_M*SPbar;    P_SA_m = q_SA_m*DPbar;       
    P_SV_M = q_SV_M*SPbar;    P_SV_m = q_SV_m*DPbar;
    P_PA_M = q_PA_M*SPbar;    P_PA_m = q_PA_m*DPbar;       
    P_PV_M = q_PV_M*SPbar;    P_PV_m = q_PV_m*DPbar;
                                                
    %% Elastances and compliances 
    
    % LV and RV minimal elastance parameters (mmHg mL^(-1))
    E_LV_m = P_LV_m / V_LV_0; 
    E_RV_m = P_RV_m / V_RV_0; 
    
    % Compliances (mL mmHg^(-1))
    C_SA = V_SA_s/P_SA_M;
    C_SV = V_SV_s/P_SV_M; 
    C_PA = V_PA_s/P_PA_M; 
    C_PV = V_PV_s/P_PV_M; 
    
    %% Resistances 

    % Arteriolar resistances (mmHg s mL^(-1)) 
    R_SA = (P_SA_M - P_SV_M)/CO;  % Systemic resistance
    R_PA = (P_PA_M - P_PV_M)/CO;  % Pulmonary resistance
    
    % Transmural resistances (mmHg s mL^(-1))
    R_tSA = 0.08; % orig 0.05
    R_tPA = 0.02; 
    
    % Valve (vlv) resistances (mmHg s mL^(-1)) 
    R_m = 2e-2; 
    R_a = 1e-4; 
    R_t = 4e-3; 
    R_p = 1e-4; 
    
    %% Heart model parameters 
    
    % Sarcomere length parameters (µm)
    Lsref   = 2; 
    Lsc0    = 1.51; 
    Lse_iso = 0.04; 
    
    % Sarcomere length shortening velocity (µm s^(-1))
    v_max   = .5*7;    
    
    % Passive stress steepness parameter  
    gamma = 7.5; % optimized from ex vivo model 
    
    %% Time-varying elastance model parameters 
    
    % Percentage of cardiac cycle 
    k_TS = 0.3; % Beginning of cardiac cycle to maximal systole  
    k_TR = 0.3; % Relaxation time fraction 
    
    %% Calculate patient-specific reference midwall surface area (Amref) for LV, SEP, and RV
    
    % Ventricular inner chamber radius (cm)
    r_LV_and_SEP = (EDV_LV * 3 / (4* pi))^(1/3); 
    r_RV         = (EDV_RV * 3 / (4* pi))^(1/3); 
    
    % Ventricle midwall radius (chamber radius (r) + 1/2 wall thickness
    % (h)) (cm)
    h_LV_and_SEP = .8; 
    h_RV         = .2; 

    % Midwall radius (cm) 
    r_m_LV_and_SEP = r_LV_and_SEP + h_LV_and_SEP/2;  
    r_m_RV         = r_RV + h_RV/2;
    
    % Outer radius (cm)
    r_o_LV_and_SEP = r_LV_and_SEP + h_LV_and_SEP; 
    r_o_RV         = r_RV + h_RV; 
    
    % Midwall reference surface area (cm^2)
    Amref_LV_and_SEP = 4 * pi * (r_m_LV_and_SEP)^2; 
    Am_RV            = 4 * pi * (r_m_RV)^2;
    
    Amref_LV  = Amref_LV_and_SEP * 2/3; % Assume LV is 2/3 of LV+SEP 
    Amref_SEP = Amref_LV_and_SEP * 1/3; % Assume SEP is 1/3 of LV+SEP
    Amref_RV  = Am_RV;
    
    %% Calculate patient-specific midwall volume (Vw) for LV, SEP, and RV 
    
    % Ventricle volume (chamber + midwall) (mL)
    Vw_chamber_LV_and_SEP = 4/3 * pi * r_o_LV_and_SEP^3;  
    Vw_chamber_RV         = 4/3 * pi * r_o_RV^3; 
    
    % Ventricular midwall volume (mL)
    Vw_LV_and_SEP = Vw_chamber_LV_and_SEP - EDV_LV; 
    Vw_RV         = Vw_chamber_RV - EDV_RV;  
    
    Vw_LV  = Vw_LV_and_SEP * 2/3; % Assume LV is 2/3 of LV+SEP 
    Vw_SEP = Vw_LV_and_SEP * 1/3; % Assume SEP is 1/3 of LV+SEP 
    
    %% Approximations for initial displacements and Amref_rv in end-diastole 
    
    % Initialize diastolic displacement values (cm)
    xm_LV_d_0  = 5; 
    xm_SEP_d_0 = 2; 
    xm_RV_d_0  = 6; 
    ym_d_0     = 3; 
    
    x0 = [xm_LV_d_0; 
        xm_SEP_d_0; 
        xm_RV_d_0;
        ym_d_0; 
        Amref_RV; 
        ]; 
    x0 = log(x0); % log-scale the initial values 
    
    % Inputs for calculating displacements 
    Vw    = [Vw_LV,Vw_SEP,Vw_RV]; 
    Amref = [Amref_LV,Amref_SEP]; 
    
    % Assume end-diastolic sarcomere length 
    SL_d    = 2; %µm 
    
    opts = optimoptions('fsolve','Display','none',...
        'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
    [fnew0,~] = fsolve(@(x) calc_xm_ym(x,Lsref,Vw,Amref,SL_d,EDV_LV,0),x0,opts); 
    fnew0 = exp(fnew0);
    
    % Outputs / Diastolic displacements
    xm_LV_d  = fnew0(1);
    xm_SEP_d = fnew0(2);
    xm_RV_d  = fnew0(3);
    ym_d     = fnew0(4);
    Amref_RV = fnew0(5); 
    
    % Set initial conditions for displacements to be used in the
    % initialconditions.m script 
    deformation.xm_LV_0  = xm_LV_d; 
    deformation.xm_SEP_0 = xm_SEP_d; 
    deformation.xm_RV_0  = xm_RV_d; 
    deformation.ym_0     = ym_d; 
    
    data.deformation = deformation; 
    
    %% Calculate passive stress parameters (k_pas) for LV and RV in end-diastole  
    
    % Midwall surface area (cm^2)
    Am_LV_d = pi * (xm_LV_d^2  + ym_d^2);
    Am_RV_d = pi * (xm_RV_d^2  + ym_d^2);
    
    % Midwall curvature (cm^(-1))
    Cm_LV_d = -2 * xm_LV_d  / (xm_LV_d^2  + ym_d^2);
    Cm_RV_d = -2 * xm_RV_d  / (xm_RV_d^2  + ym_d^2);
    
    % Midwall ratio (dimensionless) 
    z_LV_d = 3 * Cm_LV_d  * Vw_LV  / (2 * Am_LV_d); 
    z_RV_d = 3 * Cm_RV_d  * Vw_RV  / (2 * Am_RV_d); 
    
    % Instantaneous sarcomere length (µm) in end-diastole
    Ls_LV_d = SL_d; 
    Ls_RV_d = SL_d;  
    
    % Passive stress 
    sigma_pas_LV_d = (Ls_LV_d - Lsc0)^gamma; 
    sigma_pas_RV_d = (Ls_RV_d - Lsc0)^gamma; 
    
    % Dimensionless combination function
    Gamma_LV_d = -(2 / 3) * z_LV_d * (1 + (1 / 3) * z_LV_d^2 + (1 / 5) * z_LV_d^4);
    Gamma_RV_d = -(2 / 3) * z_RV_d * (1 + (1 / 3) * z_RV_d^2 + (1 / 5) * z_RV_d^4);
    
    % Passive stress scaling parameters
    k_pas_LV = EDP_LV / (Gamma_LV_d * sigma_pas_LV_d);
    k_pas_RV = EDP_RV / (Gamma_RV_d * sigma_pas_RV_d);
    
    %% Approximations for initial displacements and Amref_rv in end-systole 
    
    % Initialize systolic displacements values (cm)
    xm_LV_d_0  = 5; 
    xm_SEP_d_0 = 2; 
    xm_RV_d_0  = 6; 
    ym_d_0     = 3; 
    
    x0 = [xm_LV_d_0; 
        xm_SEP_d_0; 
        xm_RV_d_0;
        ym_d_0; 
        Amref_RV; 
        ]; 
    x0 = log(x0); % log-scale the initial values 
    
    Amref = [Amref_LV,Amref_SEP, Amref_RV]; 
    
    opts = optimoptions('fsolve','Display','none',...
        'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
    [fnew1,~] = fsolve(@(x) calc_xm_ym(x,Lsref,Vw,Amref,[],ESV_LV,1),x0,opts); 
    fnew1 = exp(fnew1);
    
    % Outputs / Systolic displacements
    xm_LV_s = fnew1(1);
    xm_RV_s = fnew1(3);
    ym_s    = fnew1(4);
    
    %% Calculate active stress parameters (k_act) for LV and RV in end-systole 
    
    % Midwall surface area (cm^2)
    Am_LV_s = pi * (xm_LV_s^2  + ym_s^2);
    Am_RV_s = pi * (xm_RV_s^2  + ym_s^2);
    
    % Midwall curvature (cm^(-1))
    Cm_LV_s = - 2 * xm_LV_s  / (xm_LV_s^2  + ym_s^2);
    Cm_RV_s = - 2 * xm_RV_s  / (xm_RV_s^2  + ym_s^2);
    
    % Midwall ratio (dimensionless)  
    z_LV_s = 3 * Cm_LV_s  * Vw_LV  / (2 * Am_LV_s); 
    z_RV_s = 3 * Cm_RV_s  * Vw_RV  / (2 * Am_RV_s); 
    
    % Myofiber strain (dimensionless)
    eps_LV_s = 0.5 * log(Am_LV_s  / Amref_LV) - (1/12) * z_LV_s^2  - 0.019 * z_LV_s^4; 
    eps_RV_s = 0.5 * log(Am_RV_s  / Amref_RV) - (1/12) * z_RV_s^2  - 0.019 * z_RV_s^4; 
    
    % Sarcomere length (µm)
    Ls_LV_s  = Lsref * exp(eps_LV_s); 
    Ls_RV_s  = Lsref * exp(eps_RV_s); 
    
    % Activation function 
    Y = .6; % set to 1 in systole 
    
    % Active stress 
    sigma_act_LV_s = Y * (Ls_LV_s  - Lsc0); 
    sigma_act_RV_s = Y * (Ls_RV_s  - Lsc0); 
    
    % Dimensionless combination function 
    Gamma_LV_s = - (2 / 3) * z_LV_s * (1 + (1 / 3) * z_LV_s^2 + (1 / 5) * z_LV_s^4);
    Gamma_RV_s = - (2 / 3) * z_RV_s * (1 + (1 / 3) * z_RV_s^2 + (1 / 5) * z_RV_s^4);
    
    % Active stress scaling parameters 
    k_act_LV = ESP_LV / (Gamma_LV_s * sigma_act_LV_s);
    k_act_RV = ESP_RV / (Gamma_RV_s * sigma_act_RV_s);


    %% Pericardium parameters
    
    Vh0 = 1.25*(EDV_LV + EDV_RV); % (mL)
    s   = 10; 
    
    %% Outputs
    
    adjpars = [C_SA; C_SV; C_PA; C_PV;                      % 1-4
        R_SA; R_tSA; R_PA; R_tPA;                           % 5-7
        R_m; R_a; R_t; R_p;                 % 8-11
        Vh0; s;                                             % 12-13
        k_pas_LV; k_pas_RV; k_act_LV; k_act_RV;             % 14-17
        ]; 
    
    pars_names = {'$C_{SA}$','$C_{SV}$','$C_{PA}$','$C_{PV}$', ...
        '$R_{SA}$','$R_{t,SA}$','$R_{PA}$','$R_{t,PA}$', ...
        '$R_{m}$','$R_{a}$','$R_{t}$','$R_{p}$',...
        '$V_{h,0}$','$s$', ...
        '$k_{pas,LV}$','$k_{pas,RV}$','$k_{act,LV}$','$k_{act,RV}$'}; 

    UB = adjpars*10; 
    LB = adjpars/10;
    
    adjpars = log(adjpars); 
    UB = log(UB);
    LB = log(LB);
    
    fixpars = [P_LV_m; P_SA_m; P_SV_m;                      % 1-3
        P_RV_m; P_PA_m; P_PV_m;                             % 4-6
        Lsref; Lsc0; Lse_iso;                               % 7-9
        v_max;                                              % 10
        V_SA_u; V_SV_u; V_PA_u; V_PV_u;                     % 11-14
        k_TS; k_TR;                                         % 15-16
        E_LV_m; E_RV_m;                                     % 17-18
        Amref_LV; Amref_SEP; Amref_RV;                      % 19-21
        Vw_LV; Vw_SEP; Vw_RV;                               % 22-24
        gamma;                                              % 25      
        ]; 
        
    data.fixpars = fixpars; 
    data.pars_names = pars_names; 

end 