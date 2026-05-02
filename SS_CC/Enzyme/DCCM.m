% Transport 3 - DCCM
% Dicarboxylate Carrier-Malate (DCCM) flux equation regulated with pH, SUC,
% Malonate, CIT, OXA
%(MAL): Pim + MALe ⇌ Pie + MALm 	
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_DCCM=DCCM(Pim,MALe,Pie,MALm,SUCm,CITm,OXAm,Hm,p)

%%% Km parameters
KA=3.33e-3;   % M Pi
KB=2.24e-3;  % M MAL
KD=KB; KC=KA;

%%% Regulation binding constant
KSUC=4.93e-04; % M SUC
KCIT=4.66e-04; % M CIT
KOXA=1.72e-03; % M OXA
KH=1.97e-06; % M H+

%Concentrations
A=Pim;
B=MALe;
C=Pie;
D=MALm;

%%% Flux
KB_prime=KB*(1+(SUCm/KSUC));
Vmax_prime = 1/(1+(Hm/KH)+(CITm/KCIT)+(OXAm/KOXA));


deno = 1+(A/KA)+(B/KB_prime)+(C/KC)+(D/KD)+(A*B/KA/KB_prime)+(C*D/KC/KD);
T_DCCM = Vmax_prime * (((A*B/KA/KB_prime))-((C*D/KC/KD)))/deno;

