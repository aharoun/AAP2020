********************************************************************************
This folder contains the replication files for the quantitative analysis

"Lack of Selection and Limits to Delegation: Firm Dynamics in Developing Countries"
by 
Ufuk Akcigit, Harun Alp and Michael Peters

OPENICPSR-119504

********************************************************************************

The folder includes the file "REPLICATION.m", that produces all results in the paper. The final output is compiled in the sub-folder "output". 

The estimation of the model is performed in the file "estimation.m". To replicate the estimation the reader has to choose which version of the model to estimate. This choice is made in the file "estimation.m" in line 10 by the variable "vers".

If vers is chosen to be empty, the baseline model is estimated:
vers = ''	 	 : Estimation of baseline model

For the robustness estimations, the following choices are available:

vers = 'WPayment'	     : Estimation using managerial payment moment (Table 9, Panel B)
vers = 'WBloomAndPayment': Estimation with managerial payment moment + Experiment (Table 9, Panel C)
vers = 'zeta04'		     : Estimation cost function elasticity = 0.4 (Table 9, Panel E)
vers = 'zeta06'		     : Estimation cost function elasticity = 0.6 (Table 9, Panel E)
vers = 'Firm'		 	 : Estimation using firm-level (instead of plant-level) data (Table 9, Panel F)
vers = 'varThetaHigh'	 : Estimation with labor supply elasticity (Table 9, Panel I)
