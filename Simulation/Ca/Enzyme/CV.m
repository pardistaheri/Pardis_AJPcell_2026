% Reaction 15: ComplexIV(CV)- Enzyme
% ComplexV flux equation
% ADPm + Pim +3Hi + Hm+  ⇌ ATPm + 3Hm+
% A - ADP; B - Pi; ; C - ATP;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

function J_CV=CV(ADPm,Pim,ATPm,dGH,pH_m,p)

%%% Thermodynamics 
dGr0=36.03; % kJ/mol Gibbs free energy of the reation at pH=7
RT=p.R_con*p.Tem;
Hm=10^(-pH_m);

betaC5= p.beta;

%%% Km parameters (Zhang 2018)
KA=0.1*50e-6; % M ADP
KB=3e-3; % M Pi
KC=50e-6; % M ATP

%%% Assign conct 
A=ADPm;
B=Pim;
C=ATPm;

%%% Flux 
deno=(1+A/KA+C/KC)*(1+B/KB);
 num=1/KA/KB*(exp(-betaC5*(dGr0-p.nH*dGH-RT*log(Hm/1e-7))/RT)*A*B-...
     exp(-(betaC5-1)*(dGr0-p.nH*dGH-RT*log(Hm/1e-7))/RT)*C);
 J_CV =1*num/deno;
