function J_FH=Flux_R(mpar,Conc)
% FUM + H2O = MAL 
% A - FUM; B - MAL;

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
drG_FH= -3.6/1000; % %standard Gibs free energy of FH refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_MAL)-(dG_FUM + dG_H2O) %dGr0 is
% clculated using the dG formation of substrates and product 
T=298.15;
pH_m=7;
Keq0=exp(-drG_FH/(R*T));
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%Binding constants 
K_A=3.5653926e-06; %FUM binding constant 
K_B=mpar(1); %MAL binding constant
Vmaxf=10*mpar(2); %Max forward reaction speed
K_H=9.0419184e-07; %proton binding constant(M) 

% metaboloite concentration in matrix
A=Conc(1); %Malm initial conc. mM
B=Conc(2); %FUMm initial conc. mM
H = 10^(-pH_m);                   

Vmaxf_prime=Vmaxf/(1+(H/K_H)+(K_H/H));

deno=(1+A/K_A)*(1+B/K_B);  
J_FH =-Vmaxf_prime/K_A/K_B*(A-(B/Keq))/deno;
end