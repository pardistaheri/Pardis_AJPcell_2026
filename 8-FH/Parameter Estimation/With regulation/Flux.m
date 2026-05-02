function J_FH=Flux(mpar,Conc)
% clculated using the dG formation of substrates and product 
R = 8.314; %gas constant (J/K/mol)
drG_FH= -3.6/1000; 
T=298;
pH_m=7.5;
Keq0=exp(-drG_FH/(R*T));
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%Binding constants 
K_A=3.5653926e-06; %FUM binding constant 
K_B=195e-06; %MALm binding constant

%Regulation binding constant
K_H=9.0419184e-07; %proton binding constant(M) 
K_ATP = mpar(1);
Vmaxf=mpar(2); %Max forward reaction speed

% metaboloite concentration in matrix
A=Conc(1); %FUMm initial conc. M
B=0; %MALm initial conc. M
H = 10^(-pH_m); 

ATP = Conc(2); %ATP conc. M

alpha =1+(ATP/K_ATP);
K_A_Prime=K_A*alpha;

Vmaxf_prime=Vmaxf/(1+(H/K_H)+(K_H/H));

deno=(1+A/K_A_Prime)*(1+B/K_B);  
J_FH =Vmaxf_prime/K_A_Prime*(A-B/Keq)/deno;
end
