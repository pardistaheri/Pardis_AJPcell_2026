% Transport 11- CU 
% Calcium Uniporter 
% Cae2+ ⇌ Cam2+
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_CU=CU(Cae,Cam,Mge,Mgm,Pie,Pim,dPsi,p);
RT=p.R_con*p.Tem;
%%% Km parameters 
% KCa=4.77e-06;
% 
% %%% Km regulations
% KMg= 2.83e-04;
% KPi= 1.91e-03;
% gamma = 4.8;
%Concentrations
A=Cae;
B=Cam;
C=Mge; %Mg_e(M)
D=Mgm; %Mg_m(M)
E=Pie; %Pi_e(M)
F=Pim; %Pi_m(M)
% 
% K_A_Prime = KCa*(1+E/KPi);K_C_Prime = KMg*(1+E/KPi);
% K_B_Prime = KCa*(1+F/KPi);K_D_Prime = KMg*(1+F/KPi);
% 
% deno = 1+(A^2/K_A_Prime^2)+(B^2/K_B_Prime^2)+(C^2/K_C_Prime^2)+(D^2/K_D_Prime^2)+((B^2*D^2)/(gamma^4*K_B_Prime^2*K_D_Prime^2))+((A^2*C^2)/(gamma^4*K_A_Prime^2*K_C_Prime^2));
% T_CU = ((A^2/K_A_Prime^2)-(B^2/K_B_Prime^2))/deno;
% 
% R = 8.314/1000; %gas constant (KJ/K/mol)
% T = 310.15; %tempreture in the isolated mito experiment(K)
FR  = 0.096484; %Faraday constant(KJ/mol/mV)
%%% Binding constants
K_A=4.7715305e-06;  
K_B=4.7715305e-06;
K_C = 5.86e-03; 
K_D = 5.86e-03;
K_E = 5.42e-04; 
K_F = 5.42e-04;

Tmax = 0.005;
nH = 2.6;
gamma = 4.8;


%%% metaboloite concentration in the experiment cell
% dPsi=140;

DeltaPhi = (2 * FR * dPsi) / (RT);
E_DeltaPhi = (DeltaPhi / nH / sinh(DeltaPhi / nH))^nH;
beta = 0.5 * (1 + log(E_DeltaPhi) / (DeltaPhi));

K_A_Prime = K_A*(1+E/K_E);K_C_Prime = K_C*(1+E/K_E);
K_B_Prime = K_B*(1+F/K_F);K_D_Prime = K_D*(1+F/K_F);

deno = 1+(A^2/K_A_Prime^2)+(B^2/K_B_Prime^2)+(C^2/K_C_Prime^2)+(D^2/K_D_Prime^2)+((B^2*D^2)/(gamma^4*K_B_Prime^2*K_D_Prime^2))+((A^2*C^2)/(gamma^4*K_A_Prime^2*K_C_Prime^2));
Tfor=((exp(beta*2*DeltaPhi)))*(A^2/K_A_Prime^2);
Trev=((exp(-(1-beta)*2*DeltaPhi)))*(B^2/K_B_Prime^2);
T_CU=Tmax*(Tfor-Trev)/deno;
