% Reaction 8: Fumarate hydratase(FH)- Enzyme
% Fumarate hydratase flux equation regulated by pH and ATP
% FUM + H2O  = MAL
% A - FUM; B - MAL;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function J_FH=FH(FUMm,MALm,ATPm,pH_m,p)
	
%%% Thermodynamics
dGr= -216.78;  % kJ/mol Gibbs free energy of the reaction at pH=7  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(pH_m-7); % pH correction 

%%% Km parameters 
KA=3.56e-06; % M FUM
KB=4.89e-05; % M MAL

%%% Regulation binding constant
KH=9.04e-7; % M H+
KATP=1.44e-5; % M ATP

%%% Concentration
A=FUMm;
B=MALm;
Hm=10^(-pH_m);

%%% flux 
alpha =1+(ATPm/KATP);
KA_prime=KA*alpha;

Vmaxf_prime=1/(1+(Hm/KH)+(KH/Hm));

deno=(1+A/KA_prime)*(1+B/KB);  
J_FH =Vmaxf_prime/KA_prime*(A-B/Keq)/deno;