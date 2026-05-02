function J_CITS=flux(mpar,Conc)
% ACOA + OXA  = COA + CIT
% A - ACOA; B - OXA; C - COA; D -CIT ;

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 298.15; %tempreture in the isolated CITS experiment(K)
pH_m = 7.4; %pH in the isolated CITS experiment(K)
drG_CITS=42.03/1000; %standard Gibs free energy of CITS refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_ACoA + dG_OXA + dG_H2O) - (dG_CoA + dG_CIT)  %dGr0 is
% clculated using the dG formation of substrates and product 
Keq0=exp(-drG_CITS/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

      
   
%%% Km parameters
K_A = 1.7784373e-05; %ACoA binding constant(M)
K_B = 1.9899582e-06; %OXA binding constant(M)
K_C = 3.3218520e-05; %CoA binding constant(M)
K_D = 3.6038156e-03; %CIT binding constant(M)

%%% Regulation binding constant
K_ATP= mpar(1);
K_ADP= mpar(2);
K_AMP= mpar(3);
K_SCOA= mpar(4);
K_H = 8.9235598e-08;

Vmaxf = 5.8469666e-01; %Maximum forward reaction rate(mmol/min)

%%% metaboloite concentration in the experiment cell
A=Conc(1); %ACoAm initial conc. M
B=Conc(2); %OXAm initial conc. M
C=0; %CoAm initial conc. M
D=0; %CITm initial conc. M
ATP=Conc(3); %ATP conc. M
ADP=Conc(4); %ADP conc. M
AMP=Conc(5); %AMP conc. M
SCOA=Conc(6); %SCOA conc. M
H_m=10^(-pH_m);

K_A_prime=K_A*(1+(ATP/K_ATP)+(ADP/K_ADP)+(AMP/K_AMP)+(SCOA/K_SCOA));
Vmaxf_prime=Vmaxf/(1+(ATP/K_ATP)+(ADP/K_ADP)+(H_m/K_H));

deno=(1+A/K_A_prime+C/K_C)*(1+B/K_B)*(1+D/K_D);  
J_CITS =Vmaxf_prime/K_A_prime/K_B*(A*B-C*D/Keq)/deno;
end