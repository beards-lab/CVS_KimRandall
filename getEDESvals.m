function [EDV, EDP, ESV, ESP] = getEDESvals(model)

    %{ 

    This function calculates the end-diastolic (ED) and end-systolic (ES)
    volumes and pressures from the model. 

    Inputs: 
    model   - output structure from model_sol.m 
    
    Ouputs: 
    EDV     - vector of end-diastolic volumes for LV and RV 
    EDP     - vector of end-diastolic pressures for LV and RV 
    ESV     - vector of end-systolic volumes for LV and RV 
    ESP     - vector of end-systolic pressures for LV and RV 

    %}

    %% Unpack outputs structure  

    % Volumes (mL)
    V_LV = model.volumes.V_LV;    
    V_RV = model.volumes.V_RV; 

    % Pressures (mmHg)
    P_LV = model.pressures.P_LV; 
    P_RV = model.pressures.P_RV; 

    %% Calculate model prediced EDP, ESP, EDV, ESV
    
    % Find max and min volumes for LV and RV 
    V_LV_min = min(round(V_LV));    V_LV_max = max(round(V_LV)); 
    V_RV_min = min(round(V_RV));    V_RV_max = max(round(V_RV));

    % EDV: max volume               % ESV: min volume
    EDV_LV = max(V_LV);             ESV_LV = min(V_LV);
    EDV_RV = max(V_RV);             ESV_RV = min(V_RV); 

    % EDP: min pressure at max volume 
    loc = find(round(V_LV) == V_LV_max); 
    EDP_LV = min(P_LV(loc));
    
    loc = find(round(V_RV) == V_RV_max); 
    EDP_RV = min(P_RV(loc));
    
    % ESP: max pressure at min volume 
    loc = find(round(V_LV) == V_LV_min);    
    ESP_LV = max(P_LV(loc));

    loc = find(round(V_RV) == V_RV_min);
    ESP_RV = max(P_RV(loc));

    %% Outputs

    EDV = [EDV_LV; EDV_RV]; 
    EDP = [EDP_LV; EDP_RV]; 
    ESV = [ESV_LV; ESV_RV]; 
    ESP = [ESP_LV; ESP_RV]; 

end 