% Transport 14- KHE 
% Potassium-proton antiporter 
% Km + He -> Ke + Hm
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_KHE=KHE(Km,He,Ke,Hm,Cam,Mgm,p)
RT=p.R_con*p.Tem;

%%% Km parameters 
KA = 7.50e-03;
KB = 9.15e-11;
KC = KA;  
KD = KB;

%%% Km regulations
KCa=1.84e-06; 
KMg=4.64e-05;

%Concentrations
A=Km; %K_m(M)
B=He; %H_e(M)
C=Ke; %K_e(M)
D=Hm; %H_m(M)

Tmax_Prime = 1/(1+(Cam/KCa)+(Mgm/KMg));

deno = 1+A/KA+C/KC+B/KB+D/KD+(A*B/(KA*KB))+((C*D)/(KC*KD));
T_KHE = Tmax_Prime*((A*B/(KA*KB))-((C*D)/(KC*KD)))/deno;

