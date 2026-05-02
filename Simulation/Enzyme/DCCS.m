% Transport 4 - DCCS
% Dicarboxylate Carrier-Succinate (DCCS) flux equation regulated with Malonate
% Pim + SUCe ⇌ Pie + SUCm 	
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_DCCS=DCCS(Pim,SUCe,Pie,SUCm,MALm,p)

%%% Km parameters
KA=3.33e-3;   % M Pi
KB=4.93e-04;  % M SUC
KD=KB; KC=KA;

%%% Regulation binding constant
KMAL=2.24e-3;  % M MAL

%Concentrations
A=Pim;
B=SUCe;
C=Pie;
D=SUCm;

%%% Flux
KB_prime=KB*(1+(MALm/KMAL));

deno = 1+(A/KA)+(B/KB_prime)+(C/KC)+(D/KD)+(A*B/KA/KB_prime)+(C*D/KC/KD);
T_DCCS = (((A*B/KA/KB_prime))-((C*D/KC/KD)))/deno;

