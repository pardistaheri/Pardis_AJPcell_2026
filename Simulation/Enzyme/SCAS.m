% Reaction 6: Succinyl COA synthetase(SCAS)- Enzyme
% Succinyl COA synthetase flux equation 
% SCOA + GDP + Pi = Suc + GTP + COA 
% A - SCOA; B - GDP; ; C - Pi; D - SUC; E - GTP; F - COA;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

function J_SCAS=SCAS(SCOAm,GDPm,Pim,SUCm,GTPm,COAm,pH_m,p)
	
%%% Thermodynamics
dGr= 1.26;  % kJ/mol Gibbs free energy of the reaction at pH=7  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(pH_m-7); % pH correction 


%%% Km parameters 
KA=7.63e-06; % M SCOA 
KB=1.83e-6; % M GDP
KC=1.81e-3; % M Pi
KD=2.35e-4; % M SUC
KE=5.17e-07; %M GTP
KF=1.25e-05; %M COA

%%% Concentration
A=SCOAm;
B=GDPm;
C=Pim;
D=SUCm;
E=GTPm;
F=COAm;

%%% flux 
deno=(1+A/KA+D/KD+F/KF+D*F/KD/KF)*(1+B/KB+C/KC+E/KE+B*C/KB/KC);
J_SCAS =(1/(KA*KB*KC))*(A*B*C-D*E*F/Keq)/deno;

