% Reaction 13: Complex III (CIII)- Enzyme 
% Complex III flux equation regulated with pH
% 2CytCo + UQH2 + 2H -> 2CytCr + UQ + 4H(+)
% A - CytCo; B - UQH2; C - CytCr; D -UQ ;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

function J_CIII=CIII(CytCoxi,UQH2m,CytCred,UQm,dPsi,dGH,pH_m,p) 

%%% Thermodynamics
dGr0=-32.53; % kJ/mol Gibbs free energy of the reaction at pH=7
RT=p.R_con*p.Tem; % 

%%% Km parameters  
KA=3.31e-06; %M CytCo 
KB=9.25e-06; %M UQH2
KC =6.34e-06; %M CytCr 
KD = 1.17e-05; %M UQ 
betaC3= p.beta;

%%% Regulation binding constant
KH1 = 3.77e-10; % M first H+
KH2 = 2.18e-07; % M second H+

%%% Assign Conct
A=CytCoxi;
B=UQH2m;
C=CytCred;
D=UQm;
Hm=10^(-pH_m);

Vmaxf_prime = 1/(KH1/Hm+1+Hm/KH2);

deno=(1+A^2/KA^2+C^2/KC^2)*(1+B/KB+D/KD);
num=1/KA^2/(KB)*(exp(-betaC3*(dGr0+4*dGH+2*RT*log(Hm/1e-7)-2*p.F_con*dPsi)/RT)*A^2*B-...
    exp(-(betaC3-1)*(dGr0+4*dGH+2*RT*log(Hm/1e-7)-2*p.F_con*dPsi)/RT)*C^2*D);
J_CIII=Vmaxf_prime*num/deno; 

