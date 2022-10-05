function venttable = makeventtable(outputs,data)

    %{ 

    This function creates a table with all of the important metrics for 
    both ventricles calculated from the model. 

    Inputs: 
    outputs     - output structure from model_sol.m 
    data        - input data structure with data and global parameters 

    Outputs: 
    venttable   - table with calculated metrics for both ventricles. 

    %}

    %% Unpack data structure 

    HR = data.HR; 

    %% Unpack outputs structure 
    
    time = outputs.time; 
    
    % Volumes (mL)
    V_LV = outputs.volumes.V_LV;    
    V_RV = outputs.volumes.V_RV; 
    
    % Pressures (mmHg)
    P_LV = outputs.pressures.P_LV;  
    P_RV = outputs.pressures.P_RV; 
    
    % Flows (L min^(-1))
    Q_a_valve = outputs.flows.Q_a_valve;    
    Q_p_valve = outputs.flows.Q_p_valve; 

    %% Calculate ventricular metrics from model outputs 

    % Stroke volume (mL) 
    SV_LV = max(V_LV) - min(V_LV); 
    SV_RV = max(V_RV) - min(V_RV); 
    
    % Ejection fraction (dimensionless) 
    EF_LV = SV_LV / max(V_LV); 
    EF_RV = SV_RV / max(V_RV);
    
    % Cardiac output (L min^(-1)) 
    CO_LV = trapz(time/60,Q_a_valve)/(time(end)/60 - time(1)/60); %SV * HR_end * 1e-3 
    CO_RV = trapz(time/60,Q_p_valve)/(time(end)/60 - time(1)/60); %SV * HR_end * 1e-3 
    
    % Cardiac power (W)  
    CP_LV = trapz(P_LV,V_LV) / 7.5 * 1e-3 * HR/60; %mean(P_sa(beat)) / 7.5 * 1e3 * SV * 1e-6 * HR_end/60  
    CP_RV = trapz(P_RV,V_RV) / 7.5 * 1e-3 * HR/60;
    
    CP_LV = CP_LV / 2; %average over 2 beats 
    CP_RV = CP_RV / 2; %average over 2 beats
    
    %% Make output table 

    Vent = ['LV';'RV']; 
    SV = [SV_LV; SV_RV]; 
    EF = [EF_LV; EF_RV];
    CO = [CO_LV; CO_RV]; 
    CP = [CP_LV; CP_RV]; 

    [EDV, EDP, ESV, ESP] = getEDESvals(outputs); 

    venttable = table(Vent,SV,EF,CO,CP,EDP,EDV,ESP,ESV);

end 

