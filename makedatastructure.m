function [data] = makedatastructure

    %{

    This function compiles all of the data into a structure.

    Inputs: 
    data    - input data structure with data and global parameters 

    Outputs: 
    data    - reassigned data structure with data from this script included

    %} 

    %% Set heart rate and heart period 

    % Heart rate (H) and heart period (T)
    HR = 60;  
    T  = 60/HR;  

    %% Set time span for the first 20 beats

    % Set tspan to solve over the first 20 beats 
    dt    = 0.001; 
    tspan = 0:dt:20*T;

    %% Pseudodata 

    % Mean systolic and diastolic blood pressure 
    SPbar = 120; 
    DPbar = 80;
    
    % Total blood volume (mL)
    Vtot = 5000; 
    
    % Cardiac Output (mL s^(-1))
    CO_lv = Vtot / 60;

    % Stroke volume (mL) 
    SV    = CO_lv*T; 
    
    % End-diastolic and end-systolic pressures and volumes 
    EDP_LV = 5;         
    EDP_RV = 5/4;           %EDP_lv * .25;  % Assume EDP_rv is 25% of EDP_lv 
    EDV_LV = 125;      
    EDV_RV = EDV_LV;        % Assume EDV_rv is equal to EDV_lv 
    ESP_LV = 1.05*SPbar;    %125;       
    ESP_RV = 25; 
    ESV_LV = EDV_LV - SV;       
    ESV_RV = ESV_LV;        % Assume EDV_rv is equal to EDV_lv 

    % Mean pulmonary arterial pressure (mmHg)
    P_PA_m = 14; 
    
    %% Set fraction of unstressed volume in the compartments 
    
    bvd_SA = .7; 
    bvd_PA = .4;
    bvd_SV = .9;
    bvd_PV = .9;
    
    %% Make data structure
    
    % Load data into a structure and convert units to SI 
    data.SPbar  = SPbar; 
    data.DPbar  = DPbar; 
    data.HR     = HR; 
    data.T      = T; 
    data.Vtot   = Vtot; 
    data.CO     = CO_lv; 
    data.SV     = SV; 
    data.tspan  = tspan; 
    data.dt     = dt;
    data.EDP_LV = EDP_LV; 
    data.EDV_LV = EDV_LV;
    data.ESP_LV = ESP_LV; 
    data.ESV_LV = ESV_LV; 
    data.EDP_RV = EDP_RV; 
    data.EDV_RV = EDV_RV; 
    data.ESP_RV = ESP_RV; 
    data.ESV_RV = ESV_RV; 
    data.P_PA_m = P_PA_m; 
    
    data.bvd_SA = bvd_SA; 
    data.bvd_PA = bvd_PA; 
    data.bvd_SV = bvd_SV; 
    data.bvd_PV = bvd_PV; 

end 