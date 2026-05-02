% Reaction 3: Aconitase(ACON)- Enzyme
% Aconitase flux equation regulated by pH
% CIT  = ICIT
% A - CIT; B - ICIT;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function J_ACON=ACON(CITm,ICITm,pH_m,p)
	
%%% Thermodynamics
dGr= 6.65;  % kJ/mol Gibbs free energy of the reaction at pH=7  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(pH_m-7); % pH correction 

%%% Km parameters 
KA=9.48e-04; % M CIT
KB=7.41e-5; % M ICIT

%%% Regulation binding constant
KH=4.27e-7; % M H+

%%% Concentration
A=CITm;
B=ICITm;
Hm=10^(-pH_m);

%%% flux 
Vmaxf_prime=1/(1+(Hm/KH));

deno=(1+A/KA)*(1+B/KB);  
J_ACON =Vmaxf_prime/KA*(A-B/Keq)/deno;