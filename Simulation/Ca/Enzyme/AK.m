% Reaction 16: Adenylate Kinase(AK)- Enzyme
% Adenylate Kinase flux equation regulated 
% 2ADP  = ATP + AMP 
% A - ADP; B - ATP; ; C - AMP; 
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

function J_AK=AK(ADPm,ATPm,AMPm,pH_m,p)
	
% %%% Thermodynamics
% dGr= -37.08;  % kJ/mol Gibbs free energy of the reaction at pH=0  
% Keq0=exp(-dGr/(p.R_con*p.Tem));
% Keq=Keq0*10^(pH_m-7); % pH correction 
%            
% %%% Km parameters 
% KA=9.63e-05; % M ADP 
% KB=2.73e-06; % M ATP
% KC=1.47e-05; % M AMP
         
%%% Concentration
A=ADPm;
B=ATPm;
C=AMPm;

%%% flux  
J_AK =(0.4331*A^2-B*C);

