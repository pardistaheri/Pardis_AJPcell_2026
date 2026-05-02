function J_MDH=flux(mpar,Conc)
% MAL + NAD = OXA + NADH + H
% A - MAL; B - NAD; ; C - OXA; D - NADH; 

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
drG_MDH= 69.46/1000; %%standard Gibs free energy of MDH refrence reaction at T=298.15 K, I=0.1, pH=0
% % dGr0 = (dG_OXA + dG_NADH)-(dG_MAL + dG_NAD) %dGr0 is
% clculated using the dG formation of substrates and product 
T=298.15;
pH_m=8;
Keq0=exp(-drG_MDH/(R*T));
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%Binding constants 
K_A=mpar(1); %MAL binding constant 
K_B=mpar(2); %NAD binding constant
K_C = mpar(3); %OXA binding constant(M)
K_D = mpar(4); %NADH binding constant (M)
Vmaxf=1.3*mpar(5); %Max forward reaction speed

% metaboloite concentration in matrix
A=Conc(1); %MALm initial conc. mM
B=Conc(2); %NADm initial conc. mM
C=Conc(3); %OXAm initial conc. mM
D=Conc(4); %NADHm initial conc. mM 

deno=(1+A/K_A)*(1+C/K_C)*(1+B/K_B+D/K_D); 
J_MDH =Vmaxf/K_A/K_B*(A*B-C*D/Keq)/deno;
end