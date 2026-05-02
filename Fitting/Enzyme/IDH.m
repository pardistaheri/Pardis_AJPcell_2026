% Reaction 4: Isocitrate dehydrogenase(IDH)- Enzyme
% Isocitrate dehydrogenase flux equation regulated by ATP, ADP,and pH
% ICIT + NAD + H2O = AKG + NADH + CO2 + 2H+
% A - ICIT; B - NAD; ; C - AKG; D - NADH; E - CO2;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function J_IDH=IDH(ICITm,NADm,AKGm,NADHm,CO2m,ATPm,ADPm,pH_m,p)
	
%%% Thermodynamics 
dGr= -4.47;  % kJ/mol Gibbs free energy of the reaction at pH=7  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(2*(pH_m-7)); % pH correction 

%%% Km parameters 
KA=4.47e-04; % M ICIT 
KB=7.21e-1; % M NAD
KC=4.46e-4; % M AKG
KD=3.13e-6; % M NADH
KE=10e-03; % M CO2
n=3.45; %Hill coefficient accounting for ICIT cooperativity 

%%% Regulation binding constant
KH=1.9e-7; % M H+
KATP=1.08e-1; % M ATP
KADP=7.07e-4; % M ADP

%%% Concentration
A=ICITm;
B=NADm;
C=AKGm;
D=NADHm;
E=CO2m;
Hm=10^(-pH_m);

%%% flux 

KA_prime=KA*(1/1+(ATPm/KATP)+(ADPm/KADP));
Vmaxf_prime=1/(1+(Hm/KH));

deno=(1+(A/KA_prime)^n)*(1+B/KB+D/KD)*(1+C/KC)*(1+E/KE);  
J_IDH =Vmaxf_prime/(KA_prime)^n/KB*((A)^n*B-C*D*E/Keq)/deno;

