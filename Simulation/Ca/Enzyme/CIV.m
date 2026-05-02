% Reaction 14: ComplexIV(CIV)- Enzyme
% ComplexIV flux equation with pH regulation
% 2CytCr + 0.5O2 + 2Hm = 2CytCo + H2O + 2dH 
% A - CytCr; B - O2; ; C - CytCo;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

function J_CIV=CIV(CytCred,O2,CytCoxi,dPsi,dGH,Cam,pH_m,p)

%%% Thermodynamics
dGr0= -122.94;  % kJ/mol Gibbs free energy of the reaction at pH=0  
RT=p.R_con*p.Tem;
betaC4= p.beta;

%%% Km parameters 
KA=1.68e-04; % M CytCr
KB=4.32e-06; % M O2
KC=1.68e-04; % M CytCo 

%%% Regulation binding constant
KH1=4.06e-06; % M 1st H+
KH2=4.87e-08; % M 2nd H+
KCa=1e-06;

%%% Concentration
A=CytCred;
B=O2;
C=CytCoxi;
Hm=10^(-pH_m);

%%% Flux 

% Vmaxf_prime = 1/(KH1/Hm+1+Hm/KH2);
Vmaxf_prime=1/((1+(5*(Cam/(KCa+Cam))))*(KH1/Hm+1+Hm/KH2));

deno=(1+A^2/KA^2+C^2/KC^2)*(1+B^0.5/KB^0.5);
num=1/KA^2/KB^.5*(exp(-betaC4*(dGr0+2*dGH+2*p.F_con*dPsi-RT*log(Hm/1e-7))/RT)*A^2*B^.5-...
    exp(-(betaC4-1)*(dGr0+2*dGH+2*p.F_con*dPsi-RT*log(Hm/1e-7))/RT)*C^2);
J_CIV=Vmaxf_prime*num/deno;
