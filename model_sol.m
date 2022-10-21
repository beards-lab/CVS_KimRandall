function [outputs,rout,J] = model_sol(adjpars,data)

    %{ 

    This function solves the time-varying in vivo version of the model to
    steady-state and then calculates 2 steady-state beats. 

    Inputs: 
    adjpars         - vector of adjustable parameters 
    data            - input data structure with data and global parameters

    Outputs: 
    outputs         - structure with all pertinent model outputs to plot 
    rout            - residual vector 
    J               - cost functional 

    %} 

    %% Unpack data structure   
    
    % Plot intermediate figures to visualize steady-state 
    plotintfigs_on = data.plotintfigs_on;  
    
    % Times (s)
    tspan = data.tspan;  
    dt    = data.dt; 
    T     = data.T; 
    
    % ODE tolerance 
    ODE_TOL = data.gpars.ODE_TOL;
    
    % Fixed parameters 
    fixpars = data.fixpars;

    %% Adjustable parameters 
    
    % Undo log from parameters.m
    adjpars = exp(adjpars); 
    
    % Compliance 
    C_SA = adjpars(1);

    %% Fixed parameters 
    
    % Unstressed volumes (mL) 
    V_SAu = fixpars(11);
    V_SVu = fixpars(12); 
    V_PAu = fixpars(13); 
    V_PVu = fixpars(14); 
    
    %% Get initial conditions
    
    init = initialconditions(adjpars,data);

    %% Set mass matrix M for DAE 

    M = speye(length(init));
    M(1,1) = 0;
    M(2,2) = 0;
    M(3,3) = 0;
    M(4,4) = 0; 

    %% Solve model 
    % Use a while loop to allow model to converge to steady-state 
    
    ndone = 0; 
    while ndone == 0 
        
        % Set options and compute model for time span tspan 
        opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL);
        sol  = ode15s(@model_invivo,[tspan(1) tspan(end)],init,opts,adjpars,data);
        
        % Plot intermediate states if plotintfigs_on is on 
        if plotintfigs_on == 1 
            % Displacements states 1-4
            figure(111)
            clf
            plot(sol.x,sol.y(1:4,:))
            legend('1','2','3','4')
            
            % Sarcomere lengths states 5-7
            figure(112)
            clf
            plot(sol.x,sol.y(5:7,:))
            legend('5','6','7')
            
            % Volumes states 8-9
            figure(113)
            clf
            plot(sol.x,sol.y(8:9,:))
            legend('8','9')
            
            % Volumes states 10-13
            figure(114)
            clf
            plot(sol.x,sol.y(10:13,:))
            legend('10','11','12','13')
        end 
    
        if sol.x(end) ~= tspan(end) 
            % Check to see if the model solved to the end of tspan. If not,
            % set the initial conditions for the next loop at the start of 
            % the previous beat and solve for 10 more beats 

            t = sol.x(1):dt:sol.x(end); 
            beats = mod(t,T); 
            x = find(round(beats,3) == 0);
            y = find(t(x) <= t(end)); 
            tspan = tspan(1):dt:tspan(x(y(end)));
        
            sols = deval(sol,tspan);
           
            if round(x(y(end-1))) == 1 
                % If model stopped after one cardiac cycle, reset initial
                % conditions to the beginning of the second cardiac cycle.
                init  = sols(:,x(y(end))); 
                tspan = tspan(x(y(end))):dt:tspan(x(y(end))) + 10*T;
            else 
                % If the model stopped after more than one cardiac cycle, 
                % reset the initial conditions to the previous fully 
                % computed cardiac cycle. 
                init  = sols(:,x(y(end-1))); 
                tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T;
            end 
             
        else 
            % If the model has successfully solved at least 10 beats, then we can
            % assess whether the model has reached steady state 
            
            % Extract sarcomere lengths and systemic arterial pressure (P_SA)     
            sols = deval(sol,tspan);
                
            Lsc_LV  = sols(5,:);
            Lsc_SEP = sols(6,:); 
            Lsc_RV  = sols(7,:); 
            P_SA    = sols(10,:) / C_SA; 
            
            % Plot extracted quantities if plotintfigs_on is on 
            if plotintfigs_on == 1
                figure(101)
                clf
                hold on 
                plot(tspan,Lsc_LV) 
                
                figure(102)
                clf
                hold on 
                plot(tspan,Lsc_SEP) 
                
                figure(103)
                clf
                hold on 
                plot(tspan,Lsc_RV) 
                
                figure(104)
                clf % new fig with each time series
                hold on 
                plot(tspan,P_SA) 
            end
            
            % Find the last 5 beats of the simulation 
            xx = find(tspan >= tspan(end) - 5*T); 
            
            % Set a peak threshold as half of the amplitude of the last 5 beats 
            threshold_LV  = (max(Lsc_LV(xx))  - min(Lsc_LV(xx)))/2;
            threshold_SEP = (max(Lsc_SEP(xx)) - min(Lsc_SEP(xx)))/2;
            threshold_RV  = (max(Lsc_RV(xx))  - min(Lsc_RV(xx)))/2;
            
            % Determine the length of half of the cardiac cycle 
            half_per = round((T/2)/dt); 
            
            % Find peaks for the last 5 beats with conditions that the
            % amplitude must be at least half of the amplitude of the last
            % 5 peats and the minimum distance between peaks must be at
            % least half of the cardiac cycle apart 
            [pks_Lsc_LV, loc_pks_Lsc_LV]  = findpeaks(...
                Lsc_LV,'MinPeakDistance',half_per,'MinPeakProminence',threshold_LV); 
            [pks_Lsc_SEP,loc_pks_Lsc_SEP] = findpeaks(...
                Lsc_SEP,'MinPeakDistance',half_per,'MinPeakProminence',threshold_SEP); 
            [pks_Lsc_RV, loc_pks_Lsc_RV]  = findpeaks(...
                Lsc_RV,'MinPeakDistance',half_per,'MinPeakProminence',threshold_RV); 
            [pks_P_SA,   loc_pks_P_SA]    = findpeaks(...
                P_SA,'MinPeakDistance',half_per); 
            
            % Exclude the last peak, so there are 4 peaks (This is
            % important in case something unexpected happens at the end of
            % the signal) 
            pks_Lsc_LV  = pks_Lsc_LV(end-5:end-1); 
            pks_Lsc_SEP = pks_Lsc_SEP(end-5:end-1); 
            pks_Lsc_RV  = pks_Lsc_RV(end-5:end-1);
            pks_P_SA    = pks_P_SA(end-5:end-1); 
            
            % Find the locations of the peaks 
            loc_pks_Lsc_LV  = loc_pks_Lsc_LV(end-5:end-1); 
            loc_pks_Lsc_SEP = loc_pks_Lsc_SEP(end-5:end-1); 
            loc_pks_Lsc_RV  = loc_pks_Lsc_RV(end-5:end-1); 
            loc_pks_P_SA    = loc_pks_P_SA(end-5:end-1); 
            
            % Find the times where the peaks occur 
            t_pks_Lsc_LV  = tspan(loc_pks_Lsc_LV);
            t_pks_Lsc_SEP = tspan(loc_pks_Lsc_SEP);
            t_pks_Lsc_RV  = tspan(loc_pks_Lsc_RV);
            t_pks_P_SA    = tspan(loc_pks_P_SA); 
            
            % Create a linear regression through the peaks 
            pf_Lsc_LV  = polyfit(t_pks_Lsc_LV,pks_Lsc_LV,1); 
            pf_Lsc_SEP = polyfit(t_pks_Lsc_SEP,pks_Lsc_SEP,1); 
            pf_Lsc_RV  = polyfit(t_pks_Lsc_RV,pks_Lsc_RV,1); 
            pf_P_SA    = polyfit(t_pks_P_SA,pks_P_SA,1); 
            
            % Extract the slope of the regression line 
            slope_Lsc_LV  = pf_Lsc_LV(1);
            slope_Lsc_SEP = pf_Lsc_SEP(1);
            slope_Lsc_RV  = pf_Lsc_RV(1);
            slope_P_SA    = pf_P_SA(1);
            
            % Plot regression line through peaks if plotintfigs_on is on 
            if plotintfigs_on == 1
                % Draw the regression line through the peaks
                y_Lsc_LV  = polyval(pf_Lsc_LV,tspan); 
                y_Lsc_SEP = polyval(pf_Lsc_SEP,tspan); 
                y_Lsc_RV  = polyval(pf_Lsc_RV,tspan); 
                y_P_SA    = polyval(pf_P_SA,tspan); 
                
                % LV sarcomere length 
                figure(101)
                hold on 
                plot(t_pks_Lsc_LV,pks_Lsc_LV,'r*')
                plot(tspan,y_Lsc_LV,'k')
                
                % SEP sarcomere length 
                figure(102)
                hold on 
                plot(t_pks_Lsc_SEP,pks_Lsc_SEP,'r*')
                plot(tspan,y_Lsc_SEP,'k')
                
                % RV sarcomere length 
                figure(103)
                hold on 
                plot(t_pks_Lsc_RV,pks_Lsc_RV,'r*')
                plot(tspan,y_Lsc_RV,'k')
                
                % Systemic arterial pressure 
                figure(104)
                hold on 
                plot(t_pks_P_SA,pks_P_SA,'r*')
                plot(tspan,y_P_SA,'k') 
            end
            
            % If the slope is sufficiently small (i.e. the line is 
            % practically horizontal), we have reached steady state 
            slope_lim = 1e-3;
                % Stopping criteria 
                if abs(slope_P_SA) < slope_lim && abs(slope_Lsc_LV) < slope_lim && ...
                        abs(slope_Lsc_SEP) < slope_lim && abs(slope_Lsc_RV) < slope_lim
                    ndone = 1; 
                end 
                
            % If we have not reached steady-state, solve the model for 10 more
            % beats and reassess convergence 
            beats = mod(tspan,T); 
            x = find(round(beats,3) == 0);
            y = find(tspan(x) <= tspan(end));
            tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T;
            init  = sols(:,x(y(end-1))); 
        end 
    end
    
    % After determining that the model is in steady-state, solve 2 more beats 
    time = [0:dt:2*T]; 
    sol  = ode15s(@model_invivo,[time(1) time(end)],init,opts,adjpars,data);
    sols = deval(sol,time);
    sols = sols'; 
    
    %% Calculate other time-varying model quantities using auxiliary equations
    
    o = zeros(37,length(time));  
    for i = 1:length(time) 
        [~,o(:,i)] = model_invivo(time(i),sols(i,:),adjpars,data);
    end 
    
    %% Outputs 
    
    outputs.time = time; 
    
    displacements.xm_LV  = sols(:,1); 
    displacements.xm_SEP = sols(:,2); 
    displacements.xm_RV  = sols(:,3); 
    displacements.ym     = sols(:,4); 
    
    lengths.Lsc_LV  = sols(:,5);
    lengths.Lsc_SEP = sols(:,6); 
    lengths.Lsc_RV  = sols(:,7); 
    
    volumes.V_LV = sols(:,8); 
    volumes.V_RV = sols(:,9); 
    volumes.V_SA = sols(:,10) + V_SAu; 
    volumes.V_SV = sols(:,11) + V_SVu; 
    volumes.V_PA = sols(:,12) + V_PAu; 
    volumes.V_PV = sols(:,13) + V_PVu; 
    volumes.Vtot = sum(sols(end,8:13)) + V_SAu + V_SVu + V_PAu + V_PVu; 
    
    pressures.P_LV = o(1,:)'; 
    pressures.P_SA = o(2,:)'; 
    pressures.P_SV = o(3,:)'; 
    pressures.P_RV = o(4,:)'; 
    pressures.P_PA = o(5,:)'; 
    pressures.P_PV = o(6,:)'; 
    
    wallvolumes.Vm_LV  = o(7,:)'; 
    wallvolumes.Vm_SEP = o(8,:)'; 
    wallvolumes.Vm_RV  = o(9,:)'; 
    
    areas.Am_LV  = o(19,:)'; 
    areas.Am_SEP = o(11,:)'; 
    areas.Am_RV  = o(12,:)'; 
    
    curvatures.Cm_LV  = o(13,:)';
    curvatures.Cm_SEP = o(14,:)';
    curvatures.Cm_RV  = o(15,:)'; 
    
    strains.eps_LV  = o(16,:)'; 
    strains.eps_SEP = o(17,:)'; 
    strains.eps_RV  = o(18,:)'; 
    
    stresses.passive.sigma_pas_LV  = o(19,:)';
    stresses.passive.sigma_pas_SEP = o(20,:)';
    stresses.passive.sigma_pas_RV  = o(21,:)';
    
    stresses.active.sigma_act_LV  = o(22,:)';
    stresses.active.sigma_act_SEP = o(23,:)';
    stresses.active.sigma_act_RV  = o(24,:)';
    
    stresses.total.sigma_LV  = o(25,:)';
    stresses.total.sigma_SEP = o(26,:)';
    stresses.total.sigma_RV  = o(27,:)';
    
    % Convert m^3 s^(-1) to L min^(-1)
    flows.Q_m = o(28,:)' * 1e-3 * 60; 
    flows.Q_a = o(29,:)' * 1e-3 * 60; 
    flows.Q_t = o(30,:)' * 1e-3 * 60; 
    flows.Q_p = o(31,:)' * 1e-3 * 60; 
    
    flows.Q_SA = o(32,:)' * 1e-3 * 60; 
    flows.Q_PA = o(33,:)' * 1e-3 * 60; 
   
    tensions.Tm_LV  = o(34,:)';
    tensions.Tm_SEP = o(35,:)'; 
    tensions.Tm_RV  = o(36,:)'; 

    activation.Y = o(37,:)'; 
    
    outputs.volumes       = volumes; 
    outputs.pressures     = pressures; 
    outputs.displacements = displacements; 
    outputs.areas         = areas;
    outputs.wallvolumes   = wallvolumes; 
    outputs.curvatures    = curvatures; 
    outputs.strains       = strains; 
    outputs.stresses      = stresses;
    outputs.lengths       = lengths; 
    outputs.flows         = flows; 
    outputs.tensions      = tensions; 
    outputs.activation    = activation; 

    % Get residual vector 
    rout = makeresidual(outputs,data); 
        
    % Calculate cost functional 
    J = rout'*rout;

end 
