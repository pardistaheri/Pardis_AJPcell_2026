function J_CII_SDH=flux(mpar,Conc)

% dGr0 = ((dfG0_Fum + dfG0_FADH2) - (dfG0_Suc + dfG0_FAD)  %dGr0 is
% clculated using the dG formation of substrates and product 
R = 8.314; %gas constant (J/K/mol)
drG_CII_SDH= -28.35/1000; 
T=298;
pH_m=7.5;
Keq0=exp(-drG_CII_SDH/(R*T));
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

%Binding constants 
K_A=1.6243574e-03; %Suc binding constant 
K_B=9.2776496e-06; %UQm binding constant
K_C=9.8995952e-02; %Fum binding constant
K_D=K_B; %UQH2 binding constant

K_OAA = mpar(1);
Vmaxf=mpar(2); %Max forward reaction speed
Vmaxr=(Vmaxf*K_C*K_D)/(Keq*K_A*K_B); %Max reverse reaction speed

K_H1 = 6.7075969e-10; %1st proton binding constant reverse reaction(M)
K_H2 = 1.6157270e-06; %2nd proton binding constant reverse reaction(M)
K_H = 1.8108068e-07; %proton binding constant forward reaction(M)

% metaboloite concentration in matrix
A=Conc(1); %Sucm initial conc. M
B=5/10^3; %UQm initial conc. M
C=0; %Fumm initial conc. M
D=0; %UQH2mm initial conc. M 
H = 10^(-pH_m);
OAA = Conc(2); %OXA conc. M

K_A_prime=K_A*(1+((OAA)/(K_OAA+OAA)));
Vmaxf_prime=Vmaxf/(1+(OAA/K_OAA)+(H/K_H));
Vmaxr_prime=Vmaxr/((K_H1/H)+1+(H/K_H2));

deno=(1+A/K_A_prime)*(1+B/K_B+D/K_D)*(1+C/K_C);  
J_CII_SDH =(Vmaxf_prime*A*B/K_A_prime/K_B-Vmaxr_prime*C*D/K_C/K_D)/deno;
end
