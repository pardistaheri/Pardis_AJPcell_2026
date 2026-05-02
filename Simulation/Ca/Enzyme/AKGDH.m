% Reaction 5: alpha-Ketoglutarate Dehydrogenase(AKGDH)- Enzyme
% alpha-Ketoglutarate Dehydrogenase flux equation regulated with ATP and
% ADP
% AKG + COA + NAD + H2O = SCOA + NADH + CO2 
% A - AKG; B - COA; ; C - NAD; D - SCOA; E - NADH; F - CO2;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

function J_AKGDH=AKGDH(AKGm,COAm,NADm,SCOAm,NADHm,CO2m,Cam,ATPm,ADPm,pH_m,p)
	
%%% Thermodynamics
dGr= -37.08;  % kJ/mol Gibbs free energy of the reaction at pH=0  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(pH_m-7); % pH correction 
           
%%% Km parameters 
KA=9.63e-05; % M AKG 
KB=2.73e-06; % M COA
KC=1.47e-05; % M NAD
KD=7.32e-06; % M SCOA
KE=4.45e-05; %M NADH
KF=10e-03; %M CO2
         
%%% Regulation binding constant
% KATP=2.23e-03; % M ATP
% KADP=3.87e-04; % M ADP
% aATP=6.84e-01; % activation ATP
% aADP=1.06e+03; % activation ADP
% KCa = 1.39e-07; %Ca binding constant
% n=2.33e-01;
% aCa=1.24e-02;

KADP = 3.90e-04; %ADP binding constant
KATP = 1.46e-03; %ATP binding constant
aADP=6.65e-01; %ADP activation factor
aATP=2.36 ; %ATP inhibition factor
% aCa=2.65; %Ca activation factor
KCa = 0.1e-06; %Ca binding constant 
aCa=20; %Ca activation factor
%%% Concentration
A=AKGm;
B=COAm;
C=NADm;
D=SCOAm;
E=NADHm;
F=CO2m;

Vmaxf = (1+((aADP*ADPm)/(ADPm+KADP)))/...
    ((1+(aATP*ATPm)/(ATPm+KATP)));

KA_prime=KA/(1+(ADPm/KADP));
KA_prime=KA_prime/((1+(Cam/KCa)+(Cam/KCa)^2)/(1+(Cam/(aCa*KCa))+(Cam/(aCa*KCa))^2));

deno=(1+(A/KA_prime))*(1+B/KB+D/KD)*(1+C/KC+E/KE)*(1+F/KF);  
J_AKGDH =(Vmaxf/(KA_prime*KB*KC))*(A*B*C-D*E*F/Keq)/deno;

