% Transport 10- Hleak 
% Passive Proton Leak 
% He+ ⇌ Hm+
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_HLEAK=HLEAK(He,Hm,dPsi,p)

RT=p.R_con*p.Tem;
%%% Km parameters (Zhang 2018)
KH=1e-7;

%Concentrations
A=He;
B=Hm;

if dPsi>1e-9
    % beta_leak is .5 for HLeak figure and 1 for the rest of the figures 
 beta_leak =1; % beta leak: OCR Sensitivity to dPsi increases with Beta
 T_HLEAK=1/KH*(A*exp(beta_leak*dPsi.*p.F_con./RT)-B*exp((beta_leak-1)*dPsi.*p.F_con./RT));
else    % no membrane potential
 T_HLEAK=0;
end

T_HLEAK =T_HLEAK;