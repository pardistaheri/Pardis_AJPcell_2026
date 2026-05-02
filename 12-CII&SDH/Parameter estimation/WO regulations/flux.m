function J_CII_SDH=flux(mpar,Conc)

% dGr0 = ((dfG0_Fum + dfG0_UQH2) - (dfG0_Suc + dfG0_UQ)  %dGr0 is
% clculated using the dG formation of substrates and product 
R = 8.314; %gas constant (J/K/mol)
drG_CII_SDH= -28.35/1000; 
T=298;
pH_m=7.5;
Keq0=exp(-drG_CII_SDH/(R*T));
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%Binding constants 
K_A=mpar(1); %Suc binding constant 
K_B=mpar(2); %UQm binding constant
K_C = 1.80*10^-3; %Fum binding constant(M)
K_D=K_B; %UQH2 binding constant

Vmaxf=mpar(3); %Max forward reaction speed
Vmaxr=(Vmaxf*K_C*K_D)/(Keq*K_A*K_B); %Max reverse reaction speed

K_H1 = 6.7075969e-10; %1st proton binding constant reverse reaction(M)
K_H2 = 1.6157270e-06; %2nd proton binding constant reverse reaction(M)
K_H = 1.8108068e-07; %proton binding constant forward reaction(M)

% metaboloite concentration in matrix
A=Conc(1); %Sucm conc. mM
B=Conc(2); %UQm conc. mM
C=Conc(3); %FUMm conc. mM
D=Conc(4); %UQH2m conc. mM 
H = 10^(-pH_m);

Vmaxf_prime=Vmaxf/(1+(H/K_H));
Vmaxr_prime=Vmaxr/((K_H1/H)+1+(H/K_H2));

deno=(1+A/K_A)*(1+B/K_B+D/K_D)*(1+C/K_C);  
J_CII_SDH =(Vmaxf_prime*A*B/K_A/K_B-Vmaxr_prime*C*D/K_C/K_D)/deno;
end
