% Transport 7- GAE
% ASP-GLU Exchanger with Ca regulation
% ASPm + GLUe ⇌ ASPe + GLUm 
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_GAE=GAE(ASPm,GLUe,ASPe,GLUm,Cam,dPsi,p)

RT=p.R_con*p.Tem;
F=p.F_con;
beta_GAE= 0.5; % ASP-GLU transporter free energy barrier

%%% Km parameters
KA=5.45e-3; % M Aspartate
KB=7.83e-3; % M Glutamate
KCa=2.25e-5; % M H+ binding constant 
KC=KA;  KD=KB;

%%%% Concnetrations
A= ASPm;
B= GLUe;
C= ASPe;
D= GLUm;

%%% Flux

Vmaxf_prime=1/(1+Cam/KCa);
deno = 1+(A/KA)+(B/KB)+(C/KC)+(D/KD)+(A*B/KA/KB)+(C*D/KC/KD);
T_GAE = Vmaxf_prime * ((exp(beta_GAE*F*dPsi/RT)*(A*B/KA/KB))-((exp(-(1-beta_GAE)*F*dPsi/RT)*(C*D/KC/KD))))/deno;

