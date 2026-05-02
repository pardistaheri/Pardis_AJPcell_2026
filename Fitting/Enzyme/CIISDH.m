% Reaction 12: Complex II-Succinate Dehydrogenase (CII-SDH)
% Complex II flux equation regulated with OXA, malonate, and pH
% SUC + UQ = FUM + UQH2
% A - SUC; B - GDP; ; C - Pi; D - SUC; E - GTP; F - COA;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

function J_CIISDH=CIISDH(SUCm,UQm,FUMm,UQH2m,OXAm,pH_m,p) 

%%% Thermodynamics 
dGr0=-2.41; % kJ/mol Gibbs free energy oreaction at pH=7
RT=p.R_con*p.Tem;
Keq0=exp(-dGr0/(RT));
Keq=Keq0; % pH correction 

% %% Km parameters 
% KA=1.62e-3; % M Succinate
% KB=9.28e-6; % M UQ 
% KC=9.9e-2; % M Fumarate
% KD=9.28e-6; % M UQH2 

%%% Assign Km parameters 
KA=1800e-6; % M Succinate
KB=140e-6; % M UQ 
KC=1800e-6; % M Fumarate
KD=2.45e-6; % M UQH2 


%%% Regulation binding constant
KH=1.81e-7; % M 1st H+ forward flux
KH1=6.71e-10; % M 1st H+ reverse flux
KH2=1.62e-6; % M 2nd H+ reverse flux
KOXA=2e-5; % M OXA

%%% Concentration
A=SUCm;
B=UQm;
C=FUMm;
D=UQH2m;
Hm=10^(-pH_m);

%%% Flux 
KA_prime=KA*(1+(OXAm/KOXA)); 

Vmaxr=(KC*KD)/(Keq*KA*KB); %Max reverse reaction speed
Vmaxf_prime=1/(1+(OXAm/KOXA)+(Hm/KH));
Vmaxr_prime=Vmaxr/((KH1/Hm)+1+(Hm/KH2));

deno=(1+A/KA_prime)*(1+B/KB+D/KD)*(1+C/KC);  
J_CIISDH =(Vmaxf_prime*A*B/KA_prime/KB-Vmaxr_prime*C*D/KC/KD)/deno;

