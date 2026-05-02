function J_GOT=flux(mpar,Conc)
% ASP + AKG =  GLU + OXA
% A - ASP; B - AKG; ; C - GLU; D - OXA;

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 303; %tempreture in the isolated PDH experiment(K)
pH_m = 7; %pH in the isolated PDH experiment(K)
drG_GOT=-1/1000; %standard Gibs free energy of SCAS refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_SUC + dG_GTP + dG_CO2)-(dG_SCOA + dG_GDP + dG_Pi) %dGr0 is
% clculated using the dG formation of substrates and product 
Keq0=exp(-drG_GOT/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant
Hm = 10^(-pH_m);

K_A = mpar(1); %ASP binding constant(M)
K_B = mpar(2); %AKG binding constant(M)
K_C = mpar(3); %GLU binding constant(M)
K_D = mpar(4); %OXA binding constant (M)
Vmaxf = mpar(5); %Maximum forward reaction rate(mmol/min)

%%% metaboloite concentration in the experiment cell
A = Conc(1); %C_ASPm (M)
B = Conc(2); %C_AKGm (M)
C = Conc(3); %C_GLUm (M) 
D = Conc(4); %C_OXAm (M)

deno=(1+A/K_A)*(1+B/K_B)*(1+C/K_C)*(1+D/K_D);
J_GOT =(Vmaxf/(K_A*K_B))*(A*B-C*D/Keq)/deno;
end