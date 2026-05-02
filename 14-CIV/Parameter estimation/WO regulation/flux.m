function J_CIV=flux(mpar,dPsi,Conc,pHi,dpH)
% 2CytCr + 0.5O2 + 2H -> 2CytCo + 2H2O + 2H(+)
% A - CytCr; B - O2; C - CytCo ;
%%% Apparent equilibrium constant
R = 8.314/1000; %gas constant (KJ/K/mol)
T = 310.15; %tempreture in the isolated mito experiment(K)
F  = 0.096484; %Faraday constant(KJ/mol/mV)
RT = R*T;
pHm = pHi+dpH; %Matrix pH in the isolated mito experiment(K)
Hi = 10^(-pHi);
Hm = 10^(-pHm);

drG_CIV=-202.2; %standard Gibs free energy of CIV refrence reaction at T=298.15 K, I=0.1, pH=0 Wu et al. 2007

%%% Binding constants for 5 parameters
K_A=mpar(1); %CytCr binding constant 
K_B=mpar(2); %O2 binding constant
K_C =K_A; %CytCo binding constant(M)
Vmaxf=mpar(3); %Max forward reaction speed (mmol/min)  

betaC4=0.3; %CIV free energy barrier
Total_C=2e-6;%total concentration of cytocrome
fract=Conc(1);

%%% metaboloite concentration in the experiment cell
A=fract*Total_C; %CytCr initial conc. mM
B=Conc(2); %O2 initial conc. mM
C=(1-fract)*Total_C; %CytCo initial conc. mM

DeltaGH=F*dPsi+RT*log(Hi/Hm);
deno=(1+A^2/K_A^2+C^2/K_C^2)*(1+B^0.5/K_B^0.5);
dGfor=exp(betaC4*((-2*DeltaGH/RT)-(2*F*dPsi/RT)-(drG_CIV/RT)));
dGrev=exp((betaC4-1)*((-2*DeltaGH/RT)-(2*F*dPsi/RT)-(drG_CIV/RT)));
J_CIV=Vmaxf/K_A^2/K_B^0.5*(dGfor*A^2*B^0.5-dGrev*C^2)/deno;

end
