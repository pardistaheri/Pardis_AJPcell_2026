function J_ACON=flux(mpar,Conc)
% CIT  = ICIT
% A - CIT; B - ICIT; 

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 298.15; %tempreture in the isolated CITS experiment(K)
pH_m = 7.4; %pH in the isolated CITS experiment(K)
drG_ACON=6.65; %standard Gibs free energy of CITS refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_CIT) - (dG_ICIT)  %dGr0 is
% clculated using the dG formation of substrates and product 
Keq0=exp(-drG_ACON/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%%% Binding constants for 5 parameters
K_A=mpar(1); %ACoA binding constant 
K_B=mpar(2); %OXA binding constant
Vmaxf=mpar(3); %Max forward reaction speed 
K_H= 4.2742102e-07;
%%% metaboloite concentration in the experiment cell
A=Conc(1); %ACoAm initial conc. mM
B=Conc(2); %OXAm initial conc. mM
H_m=10^(-pH_m);
Vmaxf_prime=Vmaxf/(1+H_m/K_H);

deno=(1+A/K_A)*(1+B/K_B);  
J_ACON =Vmaxf_prime/K_A*(A-B/Keq)/deno;
end