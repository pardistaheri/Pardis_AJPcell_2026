% Reaction 1: Pyruvate dehydrogenase (PDH)- Enzyme
% Pyruvate Dehydrogenase flux equation regulated by phosphatase and kinase
% enzymes
% PYR + COA + NAD + H2O = ACOA + CO2 + NADH + H+
% A - PYR; B - COA; C - NAD; D - ACOA; E - CO2; F - NADH;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function J_PDH=PDH(PYRm,COAm,NADm,ACOAm,CO2,NADHm,Mgm,Cam,ATPm,ADPm,pH_m,p)
	
%%% Thermodynamics
dGr= -38.75;  % kJ/mol Gibbs free energy of the reaction at pH=7  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(pH_m-7); % pH correction 

%%% Km parameters 
KA=2.7e-5; % M Pyruavte 
KB=1.21e-5; % M Coenzyme A
KC=3.52e-5; % M NAD
KD=1.21e-5; % M Acetyl-coenzyme A 
KE=10e-03; %M CO2 
KF=3.52e-5; % M NADH

%%% PDK binding constant
KHK=1.28e-7; % M H+
KATP=5.36e-5; % M ATP
KADP=1.45e-3; % M ADP
KiA=2.90e-1;  % M Pyruvate

%%% PDP binding constant
KMg=1.02e-03; % M Mg2+
KHP=4.80e-08 ; % M H+
KCa=6.280e-05; % M Ca2+
% KCa=0.1e-06; % M Ca2+
bCa=2.6; % Ca2+ activation constant

%%% Concentration
A=PYRm;
B=COAm;
C=NADm;
D=ACOAm;
E=CO2; % CO2 is a constant (No ODE for CO2) % The effect of CO2 is reflected in estmated Vmax 
F=NADHm;
Hm=10^(-pH_m);

%%% flux 
PDK = (ATPm./KATP)./((Hm./KHK+KHK./Hm).*(1+(A./KiA)).*(1 +(ATPm./KATP)+(ADPm./KADP)));
PDP = ((Mgm./KMg) + (bCa.*Cam./KCa))./((Hm./KHP+KHP./Hm).*(1 + (Mgm./KMg) + (Cam./KCa)));
Vmaxf_prime = PDP/(1+PDK);

deno=(1+A/KA)*(1+B/KB+D/KD)*(1+C/KC+F/KF)*(1+E/KE);  
J_PDH = Vmaxf_prime/KA/KB/KC*(A*B*C-D*E*F/Keq)/deno;