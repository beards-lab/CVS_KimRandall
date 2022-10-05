# CVS_KimRandall


README 

This model computes the steady-state behavior of the Kim et al model (in preparation) encompassing a closed-loop cardiovascular circuit including a transmural resistance in the systemic arterial compartment, a biventricular cardiac model accounting for ventricular interaction adapted from the Lumens et al model [1], and a patient-specific parameter value formulation. 

The codes are as follows: 

RUN_invivomodel.m: Main driver that calls code to make the data structure, load in the parameter values , call the model, and output clinical measurements and plots 

makedatastructure.m: Loads in psuedodata. Can (and should) be modified to load in patient data 

parameters.m: 	Assigns fixed parameters and calculates adjustable parameter values based on patient data 
calc_xm_ym.m:	Calculates initial estimates for displacements x_m_i and y_m for I = LW, RW, and SEP during systole and diastole 

model_sol.m: 	Solves the time-varying in vivo model to steady-state and then simulates 2 steady-state beats for analysis 
initialconditions.m: 	Calculates the initial conditions of the model accounting for consistent conditions for the DAE and potential to increase the circulating volume  
triseg.m" 	Solve the TriSeg model equations to obtain consistent conditions for the DAE variables x_m_i and y_m (repeat cardiac equations from model_invivo.m) 
model_invivo.m: 	Right hand side of the DAE consisting of a coupled cardiac and cardiovascular system model
makeresidual.m: 	Computes the residual vector 
makeventtable.m: 	Calculates routine clinical measurements and compiles them into a table 
plot_invivomodel.m: 	Plotting code 


[1] Lumens et al. Three-wall segment (TriSeg) model describing mechanics and hemodynamics of ventricular interaction. Ann Biomed Eng 37:11 2234-2255. 2009. 


![image](https://user-images.githubusercontent.com/51171638/194114034-55639c27-20b2-4bae-84a6-4880267a6458.png)
