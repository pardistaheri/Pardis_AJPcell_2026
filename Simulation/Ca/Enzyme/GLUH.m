% Transport 2 - GLUH
% Glutamate-Hydrogen co-transporter between e and m
% GLUe + He+ ⇌ GLUm + Hm+ 		
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_GLUH=GLUH(GLUe,He,GLUm,Hm,p)

%%%Km parameters (Zhang 2018)
KA=1.4e-3; % Glutamate binding constant 
KB=1e-7; % H+ binding constant 
KC=KA;  KD=KB;

%%% Concentration
A=GLUe;
B=He;
C=GLUm;
D=Hm;

%%% Flux
deno=1+A/KA+C/KC+B/KB+D/KD+A*B/KA/KB+C*D/KC/KD;
T_GLUH=1/KB/KA*(A*B-C*D)/deno;
