% Reaction 9: Malate Dehydrogenase (MDH)- Enzyme
% Malate Dehydrogenase flux equation regulated by ATP, ADP, and AMP
% MAL + NAD = OXA + NADH + H 
% A - MAL; B - NAD; ; C - OXA; D - NADH;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function J_MDH=MDH(MALm,NADm,OXAm,NADHm,ATPm,ADPm,AMPm,pH_m,p)
	
%%% Thermodynamics
dGr= 28.83;  % kJ/mol Gibbs free energy of the reaction at pH=7  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(pH_m-7); % pH correction 

%%% Km parameters 
KA=6.33e-04; % M MAL
KB=1.1e-04; % M NAD
KC=4.25e-06; % M OXA
KD=2.37e-06; % M NADH

%%% Regulation binding constant
KATP=6.95e-04; % M ATP
KADP=1.23e-03; % M ATP
KAMP=1.61e-03; % M ATP

%%% Concentration
A=MALm;
B=NADm;
C=OXAm;
D=NADHm;

%%% flux 
KD_prime=KD/(1+(ATPm/KATP)+(ADPm/KADP)+(AMPm/KAMP));

deno=(1+A/KA)*(1+C/KC)*(1+B/KB+D/KD_prime); 
J_MDH =1/KA/KB*(A*B-C*D/Keq)/deno;

