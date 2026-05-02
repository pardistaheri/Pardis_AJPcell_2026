function J_PDH=flux(mpar,Conc)
% PYR + COA + NAD (+H2O) = ACOA + CO2 + NADH
% A - PYR; B - COA; C - NAD; D -ACOA ; E - CO2tot; F - NADH;

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 310.15; %tempreture in the isolated PDH experiment(K)
pH_m = 7.5; %pH in the isolated PDH experiment(K)
drG_PDH=19.5/1000; %standard Gibs free energy of PDH refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_ACoA + dG_NADH) - (dG_PYR + dG_CoA + dG_NAD)  %dGr0 is
% clculated using the dG formation of substrates and product 
Keq0=exp(-drG_PDH/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%%% Binding constants for 4 parameters
K_A=mpar(1); %PYR binding constant 
K_B=mpar(2); %CoA binding constant
K_C=mpar(3); %NAD binding constant
K_D=K_B; %ACoA binding constant
K_E=1e-03;%CO2 binding constant
K_F=K_C; %NADH binding constant
Vmaxf=mpar(4); %Max forward reaction speed    

%%% metaboloite concentration in the experiment cell
A=Conc(1); %PYRm initial conc. mM
B=Conc(2); %COAm initial conc. mM
C=Conc(3); %NADm initial conc. mM
D=0; %ACOAm initial conc. 
E=12e-3; %CO2 initial conc. 
F=0; %NADHm initial conc. 

deno=(1+A/K_A)*(1+B/K_B+D/K_D)*(1+C/K_C+F/K_F)*(1+E/K_E);  
J_PDH =Vmaxf/K_A/K_B/K_C*(A*B*C-D*E*F/Keq)/deno;
end
