% Transport 8 - ANT (Adenine Nucleotide Translocase)
% ATP ADP Anti-transporters
% ADPe + ATPm ⇌ ADPm + ATPe 
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_ANT=ANT(ADPe,ATPm,ADPm,ATPe,dPsi,p)

RT=p.R_con*p.Tem;
F=p.F_con;
%%% Assign Km parameters 
KA=0.41*10e-6; % M ADP binding constant % modified E metelkin et al
KB=5.7*10e-6; % M ATP binding constant % modified based on E metelkin et al
KC=KA; KD=KB;


%%% Assign conct
A=ADPe;
B=ATPm;
C=ADPm;
D=ATPe;

beta_ANT=0.6; %Sensitivity increases with Beta

%%% Flux
% Shima 
deno=(1+C/KC+B/KB*exp(beta_ANT*p.F_con*dPsi/RT))*(1+A/KA+D/KD*exp((beta_ANT-1)*p.F_con*dPsi/RT));
T_ANT=((exp(beta_ANT*p.F_con*dPsi/RT)*A*B/KA/KB)-(exp((beta_ANT-1)*p.F_con*dPsi/RT)*C*D/KC/KD))/deno;

% RT=p.R_con*p.Tem;
% F=p.F_con;
% beta_ANT=0.65; %Sensitivity increases with Beta
% 
% %%% Km parameters 
% KA=4.67e-5;   % M ADPe
% KB=6.71e-03;  % M ATPm
% KC=1.08;   % M ADPm
% KD=5.70e-05;  % M ATPe
% 
% %%% Concnetration
% A=ADPe;
% B=ATPm;
% C=ADPm;
% D=ATPe;
% 
% 
% %%% Flux
% deno=(1+C/KC+((B/KB)*(exp(beta_ANT*F*dPsi/RT))))*(1+A/KA+((D/KD)*(exp((beta_ANT-1)*F*dPsi/RT))));
% Tfor=(exp(beta_ANT*F*dPsi/RT))*(A*B/KA/KB);
% Trev=(exp((beta_ANT-1)*F*dPsi/RT))*(C*D/KC/KD);
% T_ANT=(Tfor-Trev)/deno;