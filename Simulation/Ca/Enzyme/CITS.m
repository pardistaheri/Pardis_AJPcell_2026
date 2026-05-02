% Reaction 2: Citrate synthase(CTIS)- Enzyme
% Citrate Synthase flux equation regulated by ATP, ADP,AMP,SCOA,pH
% ACOA + OXA  = COA + CIT + 2H+
% A - ACOA; B - OXA; C - COA; D - CIT;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function J_CITS=CITS(ACOAm,OXAm,COAm,CITm,ATPm,ADPm,AMPm,SCOAm,pH_m,p)
	
%%% Thermodynamics
dGr= -36.5;  % kJ/mol Gibbs free energy of the reaction at pH=7  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(2*(pH_m-7)); % pH correction 

%%% Km parameters 
KA=4.99e-06; % M ACOA 
KB=4.54e-6; % M OXA
KC=4.50e-5; % M COA
KD=2.64e-3; % M CIT

%%% Regulation binding constant
KH=8.92e-8; % M H+
KATP=1.69e-3; % M ATP
KADP=2.43e-3; % M ADP
KAMP=1.31e-03; %M AMP
KSCOA=3.59e-4;  % M SCOA

%%% Concentration
A=ACOAm;
B=OXAm;
C=COAm;
D=CITm;
Hm=10^(-pH_m);

%%% flux 
KA_prime=KA*(1+(ATPm/KATP)+(ADPm/KADP)+(AMPm/KAMP)+(SCOAm/KSCOA));
Vmaxf_prime=1/(1+(ATPm/KATP)+(ADPm/KADP)+(Hm/KH));

deno=(1+A/KA_prime+C/KC)*(1+B/KB)*(1+D/KD);  
J_CITS =Vmaxf_prime/KA_prime/KB*(A*B-C*D/Keq)/deno;