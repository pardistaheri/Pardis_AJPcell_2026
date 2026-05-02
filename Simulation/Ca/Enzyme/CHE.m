% Transport 13- CHE 
% Calcium-proton antiporter 
% Cam + He -> Cae + Hm
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_CHE=CHE(Cam,He,Cae,Hm,p)
RT=p.R_con*p.Tem;

%%% Km parameters 
KA=1.72e-05;  
KB=1e-07 ;
KC = 7.95e-05; 
KD = KB;
Tmax = 0;
%Concentrations
A=Cam; %Ca_m(M)
B=He; %H_e(M)
C=Cae; %Ca_e(M)
D=Hm; %H_m(M)

deno = (1+A/KA+C/KC)+(1+(B/KB)^2+(D/KD)^2);
T_CHE =Tmax*((A*(B^2)/(KA*(KB^2)))-((C*(D^2))/(KC*(KD^2))))/deno;
