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
K_A=4.6532348e-04; %Suc binding constant 
K_B=0.0015; %FADm binding constant
K_C=0.0997; %Fum binding constant
K_D=K_B; %FADH2 binding constant
Vmaxf =mpar(3); %Max forward reaction speed
Vmaxr=(Vmaxf*K_C*K_D)/(Keq*K_A*K_B); %Max reverse reaction speed

K_H1 = 6.7075969e-10; %1st proton binding constant reverse reaction(M)
K_H2 = 1.6157270e-06; %2nd proton binding constant reverse reaction(M)
K_H = 1.8108068e-07; %proton binding constant forward reaction(M)

K_MA = mpar(1);
n=mpar(2);
% metaboloite concentration in matrix
A=Conc(1); %Sucm initial conc. M
B=8.1/10^3; %UQm initial conc. M
C=0; %Fumm initial conc. M
D=0; %UQH2m initial conc. M 
H = 10^(-pH_m);
MA = Conc(2); %MAL conc. M

K_A_prime=K_A*(1+((n*MA)/(K_MA+MA)));

Vmaxf_prime=Vmaxf/(1+(H/K_H));
Vmaxr_prime=Vmaxr/((K_H1/H)+1+(H/K_H2));

deno=(1+A/K_A_prime)*(1+B/K_B+D/K_D)*(1+C/K_C);  
J_CII_SDH =(Vmaxf_prime*A*B/K_A_prime/K_B-Vmaxr_prime*C*D/K_C/K_D)/deno;
end
