% Transport 15- NCE 
% Calcium-Sodium antiporter 
% Cam + Nae -> Cae + Nam
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_NCE=NCE(Cam,Nae,Cae,Nam,dPsi,Hm,He,p)
RT=p.R_con*p.Tem;

%%% Km parameters 
KA = 2.24e-02;
KB = 3.76e-03;
KC = 2.40e-06;  
KD = 4.99e-02;

%%% Km regulations
KHi1=2.28e-08; %proton binding constant(M)   
KHi2=6.38e-07; %proton binding constant(M)  

KH1 = 8.21e-08;
KH2 = 8.79e-08;

%Concentrations
A=Cam; %Ca_m(M)
B=Nae; %Na_e(M)
C=Cae; %Ca_e(M)
D=Nam; %Na_m(M)

KC = KC*((He/KH1)+(He^2/KH1^2)+(He^3/KH1^3)+1+(KH2/He)+(KH2^2/He^2)+(KH2^3/He^3));
KA = KA*((Hm/KH1)+1+(KH2/Hm));
Tmax_Prime = 1/((He/KHi1)+1+(KHi2/He));

deno = (1+A/KA+C/KC)+(1+(B/KB)^3+(D/KD)^3);
T_NCE = 0*((exp(0.5*p.F_con*dPsi/RT)*(A*(B^3)/(KA*(KB^3))))-(exp(-0.5*p.F_con*dPsi/RT)*((C*(D^3))/(KC*(KD^3)))))/deno;

