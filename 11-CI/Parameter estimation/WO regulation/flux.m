function J_CI=flux(mpar,Conc)
% NADH + UQ + H -> NAD + UQH2 + 4H(+)
% A - NADH; B - UQ; C - NAD; D -UQH2 ;

%%% Apparent equilibrium constant
R = 8.314/1000; %gas constant (KJ/K/mol)
T = 310.15; %tempreture in the isolated mito experiment(K)
F  = 0.096484; %Faraday constant(KJ/mol/mV)
RT = R*T;
pH_m = 8.0; %pH in the isolated CI experiment(K)
drG_CI=-109.7/1000; %standard Gibs free energy of CI refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_NADH + dG_UQ ) - (dG_NAD + dG_UQH2)  %dGr0 is clculated using the dG formation of substrates and product 
Hi = 10^(-pH_m);
Hm = 10^(-pH_m);

%%% Binding constants for 5 parameters
K_A=mpar(1); %NADH binding constant 
K_B=mpar(2); %UQ binding constant
K_C =mpar(3); %NAD binding constant(M)
K_D = mpar(4); %UQH2 binding constant(M)
K_H1=5.8057929e-07; %1st proton binding constant(M)
K_H2=1.2952258e-09; %2nd proton binding constant(M) 

Vmaxf=mpar(5); %Max forward reaction speed (mmol/min)  
betaC1=0.3; %CI free energy barrier
dPsi = 0; %mitochondrial membrane potential


%%% metaboloite concentration in the experiment cell
A=Conc(1); %NADH initial conc. mM
B=Conc(2); %UQ initial conc. mM
C=Conc(3); %NAD initial conc. mM
D=Conc(4); %UQH2 initial conc. 

DeltaGH=F*dPsi+RT*log(Hi/Hm);

Vmaxf_prime = (Vmaxf)/((Hm/K_H1)+1+(K_H2/Hm));

dGfor=exp(-betaC1*((4*DeltaGH/RT)+(drG_CI/RT)));
dGrev=exp(-(betaC1-1)*((4*DeltaGH/RT)+(drG_CI/RT)));

deno=(1+A/K_A+C/K_C)*(1+B/K_B+D/K_D);
J_CI=Vmaxf_prime/K_A/K_B*(dGfor*A*B-dGrev*C*D)/deno;
end
