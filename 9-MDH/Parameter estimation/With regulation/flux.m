function J_MDH=flux(mpar,Conc)

drG_MDH= 69.46; 
pH_m=8;
T=298.15; %K      
R  = 8.314e-3;   %gas constant [kJ/K/mol]
RT=R*T;
Keq0=exp(-drG_MDH/RT);
Keq =  Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%Binding constants 
K_A=6.33e-04; %MAL binding constant 
K_B=1.1e-04; %NAD binding constant
K_C = 4.25e-06; %OXA binding constant(M)
K_D = 2.37e-06; %NADH binding constant (M)
Vmaxf=mpar(4); %Max forward reaction speed
K_ATP = mpar(1);
K_ADP = mpar(2);
K_AMP = mpar(3);
% metaboloite concentration in matrix
A=0; %MALm initial conc. mM
B=0; %NADm initial conc. mM
C=0.1/10^3; %OXAm initial conc. mM
D=Conc(1); %NADHm initial conc. mM 
ATP = Conc(2);
ADP = Conc(3);
AMP = Conc(4);

K_D_prime = K_D/(1+(ATP/K_ATP)+(ADP/K_ADP)+(AMP/K_AMP));

deno=(1+A/K_A)*(1+C/K_C)*(1+B/K_B+D/K_D_prime); 
J_MDH =-Vmaxf/K_A/K_B*(A*B-C*D/Keq)/deno;
end
