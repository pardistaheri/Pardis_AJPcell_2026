function J_CI=flux(mpar,Conc)
% NADH + UQ + H -> NAD + UQH2 + 4H(+)
% A - NADH; B - UQ; C - NAD; D -UQH2 ;
%%% Apparent equilibrium constant
R = 8.314/1000; %gas constant (KJ/K/mol)
T = 310.15; %tempreture in the isolated mito experiment(K)
F  = 0.096484; %Faraday constant(KJ/mol/mV)
RT = R*T;
pH_m = Conc(1); %pH in the isolated CI experiment(K)
drG_CI=-109.7/1000; %standard Gibs free energy of CI refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_NADH + dG_UQ ) - (dG_NAD + dG_UQH2)  %dGr0 is clculated using the dG formation of substrates and product 
Hi = 10^(-pH_m);
Hm = 10^(-pH_m);

                        
%%% Binding constants for 5 parameters
K_A=1.7876624e-06; %NADH binding constant 
K_B=7.0410629e-05; %UQ binding constant
K_C =4.8418306e-04; %NAD binding constant(M)
K_D = 3.8913038e-04; %UQH2 binding constant(M)
betaC1=0.3; %CI free energy barrier
dPsi = 0; %mitochondrial membrane potential
K_H1=mpar(1); %1st proton binding constant(M)
K_H2=mpar(2); %2nd proton binding constant(M)
Vmaxf=mpar(3); %Max forward reaction speed(mmol/min)   

%%% metaboloite concentration in the experiment cell
A=40/10^6; %NADH conc.(M)
B=40/10^6; %UQ conc.(M)
C=0; %NAD conc.(M)
D=0; %UQH2 conc.(M)

Vmaxf_Prime = (Vmaxf)/((Hm/K_H1)+1+(K_H2/Hm));
DeltaGH=F*dPsi+RT*log(Hi/Hm);

dGfor=exp(-betaC1*((4*DeltaGH/RT)+(drG_CI/RT)));
dGrev=exp(-(betaC1-1)*((4*DeltaGH/RT)+(drG_CI/RT)));

deno=(1+A/K_A+C/K_C)*(1+B/K_B+D/K_D);
J_CI=Vmaxf_Prime/K_A/K_B*(dGfor*A*B-dGrev*C*D)/deno;
end
