function J_CIII=flux(mpar,Conc)
% 2CytCo + UQH2 + 2H -> 2CytCr + UQ + 4H(+)
% A - CytCo; B - UQH2; C - CytCr; D -UQ ;
%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 298.15; %tempreture in the isolated CIII experiment(K)
pH_m = 8.0; 
F  = 0.096484;
drG_CIII=46.69/1000; %standard Gibs free energy of CIII refrence reaction at T=298.15 K, I=0.1, pH=0

Keq0=exp(-drG_CIII/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant

   
%%% Binding constants for 5 parameters
K_A=mpar(1); %CytCo binding constant 
K_B=mpar(2); %UQH2 binding constant
K_C =mpar(3); %CytCr binding constant(M)
K_D = mpar(4); %UQ binding constant(M)
K_H1 = 3.7687838e-10; %Proton binding constant in isolated CIII(M)
K_H2 = 2.1817644e-07; %Proton binding constant in isolated CIII(M)

Vmaxf=mpar(5); %Max forward reaction speed (mmol/min)  
betaC3=0.3; %CIII free energy barrier
dPsi = 0; %mitochondrial membrane potential

A=Conc(1); %CytCo initial conc. mM
B=Conc(2); %UQH2 initial conc. mM
C=Conc(3); %CytCr initial conc. mM
D=Conc(4); %UQ initial conc. 
Hm = 10^(-pH_m);


Vmaxf_prime = Vmaxf/(K_H1/Hm+1+Hm/K_H2);

deno=(1+(A/K_A)^2+(C/K_C)^2)*(1+B/K_B+D/K_D);
dGfor=Keq.^betaC3*exp(-2*betaC3*F*dPsi/(R*T));
dGrev=Keq.^(betaC3-1)*exp(-2*(betaC3-1)*F*dPsi/(R*T));
J_CIII=Vmaxf_prime/(K_A^2)/K_B*(dGfor*(A^2)*B-dGrev*(C^2)*D)/deno;
end
