% Transport 1 - PYRH
% Pyruvate-Hydrogen co-transporter between e and m
% PYRe + He ⇌ PYRm + Hm	
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_PYR_H=PYRH(PYRe,He,PYRm,Hm,p)
a=.01;

%%% Km parameters (Zhang 2018)
KA=a*0.24e-3;    % M PYR
KB=1e-7; % M H+
KC=KA;  KD=KB;

%%% Concentration
A=PYRe;
B=He;
C=PYRm;
D=Hm;

%%% Flux
deno=1+C/KC+A/KA+D/KD+B/KB+C*D/KC/KD+A*B/KB/KA;
T_PYR_H=1/KA/KB*(A*B-C*D)/deno;
