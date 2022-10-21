%{ 
RUN_invivomodel 

Driver to run the time-varying in vivo version of the cardiac and 
cardiovascular model. 

%}

clear all

%% Flags 

% Print output figures from model 
printoutfigs_on = 0; % 0 - off; 1 - on 

% Plot intermediate figures for debugging steady-state behavior 
plotintfigs_on = 0; % 0 - off; 1 - on

%% Load pseudodata

data = makedatastructure; 

%% Get parameters 

[adjpars,~,~,data] = parameters(data);  

data.printoutfigs_on = printoutfigs_on; 
data.plotintfigs_on  = plotintfigs_on; 

%% Global parameters

ODE_TOL = 1e-8; 
data.gpars.ODE_TOL = ODE_TOL; 

%% Solve model 

outputs = model_sol(adjpars,data); 

outputtable = makeoutputtable(outputs,data)

%% Plot

plot_invivomodel(outputs,data)


