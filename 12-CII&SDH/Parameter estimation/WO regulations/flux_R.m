function J_SDH=flux_R(mpar,Conc)

% dGr0 = ((dfG0_Fum + dfG0_FADH2) - (dfG0_Suc + dfG0_FAD)  %dGr0 is
% clculated using the dG formation of substrates and product 
R = 8.314; %gas constant (J/K/mol)
drG_CII_SDH= -28.35/1000; 
T=298;
pH_m=7.5;
Keq0=exp(-drG_CII_SDH/(R*T));
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%Binding constants 
K_A=1.6244904e-03; %Suc binding constant 
K_B=mpar(1); %UQm binding constant
K_C=mpar(2); %Fum binding constant
K_D=K_B; %UQH2 binding constant
Vmaxr=mpar(3); %Max forward reaction speed

K_H1 = 6.7075969e-10; %1st proton binding constant reverse reaction(M)
K_H2 = 1.6157270e-06; %2nd proton binding constant reverse reaction(M)
% metaboloite concentration in matrix
A=Conc(1); %Sucm initial conc. mM
B=Conc(2); %UQmm initial conc. mM
C=Conc(3); %Fumm initial conc. mM
D=Conc(4); %UQH2mm initial conc. mM 
H = 10^(-pH_m);

Vmaxr_prime=Vmaxr/((K_H1/H)+1+(H/K_H2));

deno=(1+A/K_A)*(1+B/K_B+D/K_D)*(1+C/K_C);  
J_SDH =-Vmaxr_prime/K_A/K_B*(A*B-C*D/Keq)/deno;
end
