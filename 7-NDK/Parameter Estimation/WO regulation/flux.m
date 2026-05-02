function J_NDK=flux(mpar,Conc)
% GTP + ADP = GDP + ATP 
% A - GTP; B - GDP; ; C - Pi; D - SUC; E - GTP; F - COA;

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 303; %tempreture in the isolated PDH experiment(K)
pH_m = 7.4; %pH in the isolated PDH experiment(K)
drG_NDK=0/1000; %standard Gibs free energy of PDH refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_GTP + dG_ADP)-(dG_GDP + dG_ATP) %dGr0 is
% clculated using the dG formation of substrates and product 
Keq0=exp(-drG_NDK/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant
Hm = 10^(-pH_m);

K_A = mpar(1); %GTP binding constant(M)
K_B = mpar(2); %ADP binding constant(M)
K_C = K_A; %GDP binding constant(M)
K_D = K_B; %ATP binding constant (M)
Vmaxf = mpar(3); %Maximum forward reaction rate(mmol/min)
K_H = 1.0660524e-09; %proton binding constant (M)

%%% metaboloite concentration in the experiment cell
A = Conc(1); %C_GTPm (M)
B = Conc(2); %C_ADPm (M)
C = Conc(3); %C_GDPm (M) 
D = Conc(4); %C_ATPm (M)

Vmaxf_prime=Vmaxf/(1+(Hm/K_H));

deno=(1+A/K_A+C/K_C)*(1+B/K_B+D/K_D);
J_NDK =(Vmaxf_prime/(K_A*K_B))*(A*B-C*D/Keq)/deno;
end