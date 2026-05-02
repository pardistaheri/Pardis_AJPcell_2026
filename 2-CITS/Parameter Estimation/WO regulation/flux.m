function J_CITS=flux(mpar,Conc)
% ACOA + OXA  = COA + CIT
% A - ACOA; B - OXA; C - COA; D -CIT ;

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 298.15; %tempreture in the isolated CITS experiment(K)
pH_m = 8.1; %pH in the isolated CITS experiment(K)
drG_CITS=42.03/1000; %standard Gibs free energy of CITS refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_ACoA + dG_OXA + dG_H2O) - (dG_CoA + dG_CIT)  %dGr0 is
% clculated using the dG formation of substrates and product 
Keq0=exp(-drG_CITS/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%%% Binding constants for 5 parameters
K_H = 8.9235598e-08;
K_A=mpar(1); %ACoA binding constant 
K_B=mpar(2); %OXA binding constant
K_C = mpar(3); %CoA binding constant(M)
K_D = mpar(4); %CIT binding constant(M)
Vmaxf=mpar(5); %Max forward reaction speed 

%%% metaboloite concentration in the experiment cell
A=Conc(1); %ACoAm initial conc. mM
B=Conc(2); %OXAm initial conc. mM
C=Conc(3); %CoAm initial conc. mM
D=Conc(4); %CITm initial conc. 
H_m=10^(-pH_m);
Vmaxf_prime=Vmaxf/(1+(H_m/K_H));

deno=(1+A/K_A+C/K_C)*(1+B/K_B)*(1+D/K_D);  
J_CITS =Vmaxf_prime/K_A/K_B*(A*B-C*D/Keq)/deno;
end