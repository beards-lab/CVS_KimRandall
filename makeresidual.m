function [rout] = makeresidual(outputs,data)

    %{

    This function builds the residual vector for the local sensitivity
    analysis and optimization schemes. 

    Inputs: 
    outputs     - output structure from model_sol.m 
    data        - input data structure with data and global parameters 
    
    Outputs: 
    rout        - residual vector 

    %} 

    %% Unpack data structure 
    
    EDV_LV_data  = data.EDV_LV; 
    ESV_LV_data  = data.ESV_LV; 
    EDV_RV_data  = data.EDV_RV; 
    ESV_RV_data  = data.ESV_RV;
    
    CO_data     = data.CO/1000*60; % L / min (input data) 
    DPbar_data  = data.DPbar; 
    SPbar_data  = data.SPbar; 
    P_PA_m_data = data.P_PA_m; % mmHg 
    
    %% Unpack outputs structure 
    time = outputs.time; 

    V_LV = outputs.volumes.V_LV; 
    V_RV = outputs.volumes.V_RV; 

    P_SA = outputs.pressures.P_SA; 
    P_PA = outputs.pressures.P_PA;

    Q_a_valve = outputs.flows.Q_a_valve;

    %% Calculate metrics 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cardiac output (L min^(-1))
    CO_model = trapz(time/60,Q_a_valve)/(time(end)/60 - time(1)/60); 
    rout_c1 = (CO_model - CO_data)/CO_data; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LV end-diastolic and end-systolic volumes (mL) 
    EDV_LV_model = max(V_LV); 
    ESV_LV_model = min(V_LV); 
    rout_v1 = [(EDV_LV_model - EDV_LV_data)/EDV_LV_data; 
        (ESV_LV_model - ESV_LV_data)/ESV_LV_data]; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RV end-diastolic and end-systolic volumes (mL) 
    EDV_RV_model = max(V_RV); 
    ESV_RV_model = min(V_RV); 
    rout_v2 = [(EDV_RV_model - EDV_RV_data)/EDV_RV_data; 
        (ESV_RV_model - ESV_RV_data)/ESV_RV_data];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diastolic systemic arterial pressure (mmHg)
    DPbar_model = min(P_SA); 
    rout_p1 = (DPbar_model - DPbar_data)/DPbar_data; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Systolic systemic arterial pressure (mmHg)
    SPbar_model = max(P_SA); 
    rout_p2 = (SPbar_model - SPbar_data)/SPbar_data; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean pulmonary arterial pressure (mmHg)
    P_PA_m_model = trapz(time,P_PA)/(time(end) - time(1));
    rout_p3 = (P_PA_m_model - P_PA_m_data)/P_PA_m_data;   

    %% Outputs
    
    rout = [rout_c1;  
        rout_v1; rout_v2; 
        rout_p1; rout_p2; rout_p3; 
        ]; 

end 