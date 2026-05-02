function J_SCAS=flux(mpar,Conc)
% SCOA + GDP + Pi = Suc + GTP + COA 
% A - SCOA; B - GDP; ; C - Pi; D - SUC; E - GTP; F - COA;

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 303; %tempreture in the isolated PDH experiment(K)
pH_m = 7.4; %pH in the isolated PDH experiment(K)
drG_SCAS=47.62/1000; %standard Gibs free energy of SCAS refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_SUC + dG_GTP + dG_CO2)-(dG_SCOA + dG_GDP + dG_Pi) %dGr0 is
% clculated using the dG formation of substrates and product 
Keq0=exp(-drG_SCAS/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant
Hm = 10^(-pH_m);

K_A = mpar(1); %SCOA binding constant(M)
K_B = mpar(2); %GDP binding constant(M)
K_C = mpar(3); %Pi binding constant(M)
K_D = mpar(4); %SUC binding constant (M)
K_E = mpar(5); %GTP binding constant(M)
K_F = mpar(6); %COA binding constant (M)
Vmaxf = mpar(7); %Maximum forward reaction rate(mmol/min)


%%% metaboloite concentration in the experiment cell
A = Conc(1); %C_SCOAm (M)
B = Conc(2); %C_GDPm (M)
C = Conc(3); %C_Pim (M) 
D = Conc(4); %C_SUCm (M)
E = Conc(5); %C_GTPm (M) 
F = Conc(6); %C_COAm (M)

deno=(1+A/K_A+D/K_D+F/K_F+D*F/K_D/K_F)*(1+B/K_B+C/K_C+E/K_E+B*C/K_B/K_C);
J_SCAS =(Vmaxf/(K_A*K_B*K_C))*(A*B*C-D*E*F/Keq)/deno;
end