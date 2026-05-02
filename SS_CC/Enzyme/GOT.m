% Reaction 5: Glutamate Oxaloacetate Transaminase(GOT)- Enzyme
% Glutamate Oxaloacetate Transaminase flux equation 
% ASP + AKG =  GLU + OXA
% A - ASP; B - AKG; ; C - GLU; D - OXA;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

function J_GOT=GOT(ASPm,AKGm,GLUm,OXAm,pH_m,p)
	
%%% Thermodynamics
dGr= -1.47;  % kJ/mol Gibbs free energy of the reaction at pH=7  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(pH_m-7); % pH correction 

            
%%% Km parameters 
KA=2.11e-03; % M ASP 
KB=2.19e-04; % M AKG
KC=1.45e-02; % M GLU
KD=1.07e-04; % M OXA

% %%% Assign Km parameters (Zhang 2018)
% KA=3.9e-3;  % M Asparate % xiao comment use jason's parameters
% KB=430e-6;   % M alpha-ketoglutrate 
% KC=8.9e-3;   % M Glutamate 
% KD=88e-6;    % M Oxaloacetate 


%%% Concentration
A=ASPm;
B=AKGm;
C=GLUm;
D=OXAm;

%%% flux 
deno=(1+A/KA)*(1+B/KB)*(1+C/KC)*(1+D/KD);
J_GOT =(1/(KA*KB))*(A*B-C*D/Keq)/deno;

