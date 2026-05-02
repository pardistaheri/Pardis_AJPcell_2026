% Transport 9 - PIC (Inorganic Phosphate Carrier)
% Phosphate-H cotransporter
% Pie + He+ ⇌ Pim + Hm+
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_PiH_C=PIC(Pie,He,Pim,Hm,p)

%%% Km parameters (Zhang 2018)
KA=9.4e-3; % Phosphate binding constant 
KB=1e-7; % H+ binding constant 
KC=KA; KD=KB;

%%% Assign conct
A=Pie;
B=He;
C=Pim;
D=Hm;

%%% Flux
deno=1+C/KC+A/KA+D/KD+B/KB+C*D/KC/KD+A*B/KB/KA;
T_PiH_C=1/KB/KA*(A*B-C*D)/deno;
