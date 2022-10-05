function [] = plot_invivomodel(outputs,data)


    %% Unpack outputs structure 
    time = outputs.time; 
    
    % Volumes (mL)
    V_LV = outputs.volumes.V_LV;    V_RV = outputs.volumes.V_RV; 
    
    % Pressures (mmHg)
    P_LV = outputs.pressures.P_LV;  P_RV = outputs.pressures.P_RV; 
    P_SA = outputs.pressures.P_SA;  P_PA = outputs.pressures.P_PA; 
    P_SV = outputs.pressures.P_SV;  P_PV = outputs.pressures.P_PV; 
    
    % Flows (L min^(-1))
    Q_m_valve = outputs.flows.Q_m_valve;    Q_t_valve = outputs.flows.Q_t_valve; 
    Q_a_valve = outputs.flows.Q_a_valve;    Q_p_valve = outputs.flows.Q_p_valve; 
    
    % Sarcomere lengths (um) 
    Lsc_LV  = outputs.lengths.Lsc_LV; 
    Lsc_SEP = outputs.lengths.Lsc_SEP; 
    Lsc_RV  = outputs.lengths.Lsc_RV; 

    % Septal curvature (cm^{-1}) 
    Cm_SEP = outputs.curvatures.Cm_SEP; 

    y_v = outputs.activation.y_v; 

    %% Unpack data structure 

    SPbar = data.SPbar; 
    DPbar = data.DPbar;

    printoutfigs_on = data.printoutfigs_on; 
    
    %% Make EDPVR Klotz curve 
    
    V_EDPVR = [10:300]; 
    [EDV, EDP] = getEDESvals(outputs); 
    [P_LV_EDPVR,P_RV_EDPVR] = makeKlotzcurve(EDV,EDP,V_EDPVR); 

    %% Find times for end-systole and end-diastole 

    i_ES = find(diff(Q_m_valve) > 0,1,'first'); 
    i_ED = find(diff(Q_a_valve) > 0,1,'first'); 
    
    ES = outputs.time(i_ES); 
    ED = outputs.time(i_ED); 

    %% Figures 

    purple = [148, 0, 211]/255; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PV loop 
    hfig1 = figure(1); 
    clf
    hold on 
    h1 = plot(V_LV, P_LV, 'r','linewidth',2);
    h2 = plot(V_RV, P_RV, 'b','linewidth',2);
    plot(V_EDPVR,P_LV_EDPVR,'r')
    plot(V_EDPVR,P_RV_EDPVR,'b')
    xlabel('Volume (mL)')
    ylabel('Pressure (mmHg)')
    legend([h1 h2],'LV','RV')

    x_max = max(vertcat(V_LV,V_RV));
    x_max = max(x_max,150); 
    y_max = max(vertcat(P_LV,P_RV)); 
    xlim([0,round(x_max)+10])
    ylim([0,round(y_max)+10])
    set(gca,'FontSize',20)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Psa vs time 
    hfig2 = figure(2);  
    clf
    hold on 
    plot([time(1) time(end)],SPbar * ones(2,1),'k:','linewidth',0.5)
    plot([time(1) time(end)],DPbar * ones(2,1),'k:','linewidth',0.5)
    h1 = plot(time,P_SA, 'color', purple, 'linewidth',2); 
    xlabel('Time (s)')
    ylabel('Pressure (mmHg)','linewidth',2)
    legend([h1],'P_{SA}')
    ylim([round(DPbar)-20 round(SPbar)+20])
    set(gca,'FontSize',20)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % High pressures vs time 
    hfig3 = figure(3);
    clf 
    hold on
    plot(time,P_LV,'r','linewidth',2)
    plot(time,P_SA,'color',purple,'linewidth',2)
    yline(data.SPbar,':')
    yline(data.DPbar,':')
    legend('P_{LV}','P_{SA}','orientation','horizontal')
    xlabel('Time (s)')
    ylabel('Pressure (mmHg)')
    
    y_max = max(vertcat(P_LV,P_SA)); 
    xlim([time(1) time(end)])
    ylim([0 round(y_max)+10])
    set(gca,'FontSize',20)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Low pressures vs time 
    hfig4 = figure(4);
    clf 
    hold on
    plot(time,P_SV,'k','linewidth',2)
    plot(time,P_RV,'b','linewidth',2)
    plot(time,P_PA,'c','linewidth',2)
    plot(time,P_PV,'m','linewidth',2)
    legend('P_{SV}','P_{RV}','P_{PA}','P_{PV}','orientation','horizontal')
    xlabel('Time (s)')
    ylabel('Pressure (mmHg)')
    set(gca,'FontSize',20)
    xlim([time(1) time(end)])
    ylim([0 max(P_PA + 10)])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ventricular volumes vs time 
    hfig5 = figure(5); 
    clf
    hold on 
    plot(time,V_LV,'r','linewidth',2)
    plot(time,V_RV,'b','linewidth',2)
    yline(data.EDV_LV,':')
    yline(data.ESV_LV,':')
    legend('LV','RV','orientation','horizontal')
    xlabel('Time (s)')
    ylabel('Volume (mmHg)')
    set(gca,'FontSize',20)
    xlim([time(1) time(end)])
    if max(V_LV) >= max(V_RV)
        maxV = max(V_LV);
    else 
        maxV = max(V_RV);
    end
    if min(V_LV) <= min(V_RV)
        minV = min(V_LV);
    else 
        minV = min(V_RV);
    end
    ylim([(minV-50) (maxV+50)])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Valve flows vs time 
    hfig6 = figure(6); 
    clf
    hold on 
    plot(time,Q_t_valve,'b','linewidth',2)
    plot(time,Q_p_valve,'b--','linewidth',2)
    plot(time,Q_a_valve,'r--','linewidth',2)
    plot(time,Q_m_valve,'r','linewidth',2)
    legend('Q_{t}','Q_{p}','Q_{a}','Q_{m}')
    xlabel('Time (s)')
    ylabel('Flow (L min^{-1})')
    set(gca,'FontSize',20)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sarcomere lenghts vs time 
    hfig7 = figure(7); 
    clf
    hold on 
    plot(time,Lsc_LV,'r')
    plot(time,Lsc_SEP,'g')
    plot(time,Lsc_RV,'b')
    legend('LV','SEP','RV')
    xlabel('Time (s)')
    ylabel('Sarcomere length (\mum)')
    set(gca,'FontSize',20)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Activation function 
    hfig8 = figure(8); 
    clf
    plot(time,y_v,'k','linewidth',3)
    set(gca,'FontSize',20)
    xlabel('Time (s)')
    ylabel('Y(t)')
    title('Cardiac Activation')

    xline(ES,':')
    xline(ES+data.T,':')
    xline(ED,'--')
    xline(ED+data.T,'--')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Septal curvature 
    hfig9 = figure(9); 
    clf
    hold on
    plot(time,Cm_SEP,'k','linewidth',2) 
    ylabel('Septal Curvature (cm^{-1})')
    xlabel('Time (s)')
    set(gca,'FontSize',20)
    ylim([0 .5])

    xline(ES,':')
    xline(ES+data.T,':')
    xline(ED,'--')
    xline(ED+data.T,'--')


    %% Print figures

    if printoutfigs_on == 1
        if ~exist('Figures', 'dir')
            mkdir('Figures')
        end

        print(hfig1,'-dpng','Figures/F1_PVloops.png')
        print(hfig2,'-dpng','Figures/F2_Psa.png')
        print(hfig3,'-dpng','Figures/F3_HighP.png')
        print(hfig4,'-dpng','Figures/F4_LowP.png')
        print(hfig5,'-dpng','Figures/F5_Volume.png')
        print(hfig6,'-dpng','Figures/F6_Fvlv.png')
        print(hfig7,'-dpng','Figures/F7_SL.png')
        print(hfig8,'-dpng','Figures/F8_Activation.png')
        print(hfig9,'-dpng','Figures/F9_Curvature.png')

    end 

end 
