function J_AKGDH=flux(mpar,Conc)
% AKG + COA + NAD + H2O = SCOA + NADH + CO2 + H
% A - AKG; B - COA; ; C - NAD; D - SCOA; E - NADH; F - CO2;

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 303; %tempreture in the isolated PDH experiment(K)
pH_m = 7.2; %pH in the isolated PDH experiment(K)
drG_AKGDH=12.82/1000; %standard Gibs free energy of PDH refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_AKG + dG_COA + dG_H2O)-(dG_SCOA + dG_NADH + dG_CO2) %dGr0 is
% clculated using the dG formation of substrates and product 
Keq0=exp(-drG_AKGDH/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant
Hm = 10^(-pH_m);

K_A = 1.58e-03; %AKG binding constant(M)
K_B = 2.73e-06; %COA binding constant(M)
K_C = 1.47e-05; %NAD binding constant(M)
K_D = 7.32e-06; %SCOA binding constant(M)
K_E = 4.45e-05; %NADH binding constant(M)
K_F = 1e-03; %CO2 binding constant
Vmaxf = mpar(5); %Maximum forward reaction rate(mmol/min)

K_ADP = mpar(1); %ADP binding constant
K_ATP = mpar(2); %ATP binding constant
alpha_ADP=mpar(3);
alpha_ATP=mpar(4);


%%% metaboloite concentration in the experiment cell
A = Conc(1); %C_AKGm (M)
B = Conc(2); %C_COAm (M)
C = Conc(3); %C_NADm (M) 
D = 0; %C_SCOAm (M)
E = 0; %C_NADHm (M) 
F=12e-3; %CO2m conc. M

ADP=Conc(4); %ADPm conc. 
ATP=Conc(5); %ATPm conc. 


Vmaxf_prime = Vmaxf*((1+(alpha_ADP*ADP)/(ADP+K_ADP)))/...
    ((1+(alpha_ATP*ATP)/(ATP+K_ATP)));

K_A=K_A/((1+ADP/K_ADP));

deno=(1+(A/K_A))*(1+B/K_B+D/K_D)*(1+C/K_C+E/K_E)*(1+F/K_F);  
J_AKGDH =(Vmaxf_prime/(K_A*K_B*K_C))*(A*B*C-D*E*F/Keq)/deno;
end