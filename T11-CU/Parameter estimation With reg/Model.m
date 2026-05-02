function J_CU=Model(mpar,Conc)
% Cae -> Cam
R = 8.314/1000; %gas constant (KJ/K/mol)
T = 310.15; %tempreture in the isolated mito experiment(K)
FR  = 0.096484; %Faraday constant(KJ/mol/mV)
%%% Binding constants
K_A=4.7715305e-06;  
K_B=4.7715305e-06;
K_C = mpar(1); 
K_D = mpar(1);
K_E = mpar(2); 
K_F = mpar(2);

Tmax = mpar(3);
nH = 2.6;
gamma = 4.8;


%%% metaboloite concentration in the experiment cell
A=Conc(1); %Ca_e(M)
B=0; %Ca_m(M)
C=Conc(2); %Mg_e(M)
D=0; %Mg_e(M)
E=Conc(3); %Pi_e(M)
F=0; %Pi_e(M)
dPsi=Conc(4);

DeltaPhi = (2 * FR * dPsi) / (R * T);
E_DeltaPhi = (DeltaPhi / nH / sinh(DeltaPhi / nH))^nH;
beta = 0.5 * (1 + log(E_DeltaPhi) / (DeltaPhi));

K_A_Prime = K_A*(1+E/K_E);K_C_Prime = K_C*(1+E/K_E);
K_B_Prime = K_B*(1+F/K_F);K_D_Prime = K_D*(1+F/K_F);

deno = 1+(A^2/K_A_Prime^2)+(B^2/K_B_Prime^2)+(C^2/K_C_Prime^2)+(D^2/K_D_Prime^2)+((B^2*D^2)/(gamma^4*K_B_Prime^2*K_D_Prime^2))+((A^2*C^2)/(gamma^4*K_A_Prime^2*K_C_Prime^2));
Tfor=((exp(beta*2*DeltaPhi)))*(A^2/K_A_Prime^2);
Trev=((exp(-(1-beta)*2*DeltaPhi)))*(B^2/K_B_Prime^2);
J_CU=Tmax*(Tfor-Trev)/deno;
end
