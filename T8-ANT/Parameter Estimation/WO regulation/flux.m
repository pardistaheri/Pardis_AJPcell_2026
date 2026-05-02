function J_ANT=flux(mpar,Conc)
% ADPe + ATPm -> ADPm + ATPe
% A - ADPe; B - ATPm; C - ADPe ; D - ATPm ;
%%% Apparent equilibrium constant
R = 8.314/1000; %gas constant (KJ/K/mol)
T = 310.15; %tempreture in the isolated mito experiment(K)
F  = 0.096484; %Faraday constant(KJ/mol/mV)
RT = R*T;

%%% Binding constants for 5 parameters
K_A=mpar(1); %ADPe binding constant 
K_B=mpar(2); %ATPm binding constant
K_C= mpar(3); %ADPm binding constant 
K_D= mpar(4); %ATPe binding constant
Tmaxf=mpar(5); %Max forward reaction speed (mmol/min)  

beta_ANT=0.65; %CIV free energy barrier

%%% metaboloite concentration in the experiment cell
A=Conc(1); %ADP_free conc. M outside of MM
B=Conc(2); %ATP_free conc. M inside of matrix
C=Conc(3); %ADP_free conc. M inside of matrix
D=Conc(4); %ATP_free conc. M outside of MM
dPsi=Conc(5);

deno=(1+C/K_C+((B/K_B)*(exp(beta_ANT*F*dPsi/RT))))*(1+A/K_A+((D/K_D)*(exp((beta_ANT-1)*F*dPsi/RT))));
Tfor=(exp(beta_ANT*F*dPsi/RT))*(A*B/K_A/K_B);
Trev=(exp((beta_ANT-1)*F*dPsi/RT))*(C*D/K_C/K_D);
J_ANT=Tmaxf*(Tfor-Trev)/deno;
end
