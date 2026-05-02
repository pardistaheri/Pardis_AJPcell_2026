% Reaction 11: Complex I (CI)- Enzyme  
% Complex I flux equation regulated with Ca2+ and pH
% NADH + UQ + H = NAD + UQH2 + 4dH 
% A - NADH; B - UQ; ; C - NAD; D - UQH2;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

function J_CI=CI(NADHm,UQm,NADm,UQH2m,Cam,dGH,pH_m,p)

%%% Thermodynamics
dGr0=-69.37; % kJ/mol Gibbs free energy of the reaction at pH=7
RT=p.R_con*p.Tem; % 

betaC1= p.beta;

%%% Km parameters
KA=1.79e-6; % M NADH
KB=7.01e-5; % M UQ
KC=4.93e-4; % M NAD
KD=3.96e-4; % M UQH2

%%% Regulation binding constant
KH1=5.81e-7; % M 1st H+
KH2=1.29e-9; % M 2nd H+

% KCa=3.34e-05; % M Ca2+
KCa=5e-06; % M Ca2+
%%% Concentration
A=NADHm;
B=UQm;
C=NADm;
D=UQH2m;
Hm=10^(-pH_m);

%%% Flux
% KA=KA*(1+((Cam)/(KCa+Cam)));

% Vmaxf_prime=1/(1+(Cam/(KCa+Cam))+(Hm/KH1)+(KH2/Hm));
Vmaxf_prime=1/(1+(100*(Cam/(KCa+Cam)))*((Hm/KH1)+1+(KH2/Hm)));

deno=(1+A/KA+C/KC)*(1+B/KB+D/KD);
num=1/KA/KB*(exp(-betaC1*(dGr0+4*dGH-RT*log(Hm/1e-7))/RT)*A*B-...
    exp(-(betaC1-1)*(dGr0+4*dGH-RT*log(Hm/1e-7))/RT)*C*D);
J_CI=Vmaxf_prime*num/deno;


