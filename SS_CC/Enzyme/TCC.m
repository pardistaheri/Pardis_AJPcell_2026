% Transport 5 - TCC
% MAL-HCIT  Antiporter
% Tricarboxylate Carrier (TCC)
% TCC: MALe + HCITm ⇌ MALm	+ HCITe
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_TCC=TCC(CITm,Hm,MALe,CITe,He,MALm,p)

%%%Km parameters (Zhang 2018)
KA=1e-3;   % M CIT binding constant
KB=0.25e-3;   % M Malate binding constant
KH=1e-7;  % M H+ binding constant 
KC=KA;  KD=KB;

%%% Assign conct
A=Hm*CITm; % HCITm
B=MALe;
C=He*CITe; % HCITe
D=MALm;

%%% Flux
deno=1+A/KA/KH+B/KB+C/KC/KH+D/KD+A*B/KA/KB/KH+C*D/KC/KD/KH;
T_TCC=1/KB/KA/KH*(A*B-C*D)/deno;
