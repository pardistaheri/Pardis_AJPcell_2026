% Transport 12- NHE 
% Sodium-proton antiporter 
% Nam + He -> Nae + Hm
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_NHE=NHE(Nam,He,Nae,Hm,p)
RT=p.R_con*p.Tem;
%%% Km parameters 
KA=2.61e-04;  
KB=1.003e-09 ;
KC = KA; 
KD = KB;

%%% Km regulations
KiH = 6.97e-08;
nh = 3.45;

%Concentrations
A=Nam; %Na_m(M)
B=He; %H_e(M)
C=Nae; %Na_e(M)
D=Hm; %H_m(M)

Tmax_Prime = D.^nh/(KiH.^nh+D.^nh);
deno = (1+A/KA+C/KC)+(1+B/KB+D/KD);
T_NHE = Tmax_Prime*(A*B/KA/KB-C*D/KC/KD)/deno;

