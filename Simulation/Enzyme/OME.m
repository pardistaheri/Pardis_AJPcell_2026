% Transport 6 - OME
% AKG-MAL exchanger (OME)
% AKGe + MALm ⇌ AKGm	+ MALe
% A - AKGe; B - MALm; C - AKGm; D - MALe;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_OME=OME(AKGe,MALm,AKGm,MALe,p)

%%% Km parameters
KA=2.26e-4;   % M AKG binding constant
KB=1.28e-3;   % M MAL binding constant
KC=KA;  KD=KB;

%%% Concnetrations
A=AKGe; 
B=MALm;
C=AKGm; 
D=MALe;

%%% Flux
deno = 1+(A/KA)+(B/KB)+(C/KC)+(D/KD)+(A*B/KA/KB)+(C*D/KC/KD);
T_OME =((A*B/KA/KB)-(C*D/KC/KD))/deno;
