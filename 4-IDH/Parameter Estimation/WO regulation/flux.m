function J_IDH=flux(mpar,Conc)
% ICIT+ NAD+ H2O↔ CO2+ NADH+ AKG% A - ICIT; B - NAD; C - NADH;

%%% Apparent equilibrium constant
R = 8.314; %gas constant (J/K/mol)
T = 303; %tempreture in the isolated PDH experiment(K)
pH_m = 7.2; %pH in the isolated PDH experiment(K)
drG_IDH=91.77/1000; %standard Gibs free energy of PDH refrence reaction at T=298.15 K, I=0.1, pH=0
% dGr0 = (dG_AKG + dG_NADH + dG_CO2) - (dG_ICIT + dG_NAD + dG_H2O)  %dGr0 is
% clculated using the dG formation of substrates and product 
Keq0=exp(-drG_IDH/(R*T)); %standard Keq
Keq=Keq0*10^(pH_m); %pH-dependent and binding-dependant apparent equilibrium constant
Hm = 10^(-pH_m);


K_iH = 1.9e-07; %Proton binding constant in isolated IDH(M)

K_A = 4.4715300e-04; %ICIT binding constant(M)
K_B=mpar(2); %NAD binding constant
K_C=mpar(3); %AKG binding constant
K_D=mpar(4); %NADH binding constant
K_E=1e-03;%CO2 binding constant
Vmaxf=mpar(5); %Max forward reaction speed  
n=3.45;

%%% metaboloite concentration in the experiment cell
A=Conc(1); %ICITm conc. M
B=Conc(2); %NADm conc. M
C=0; %AKGm conc. M
D=Conc(3); %NADHm conc. M
E=12e-3; %CO2m conc. M

Vmaxf_prime = Vmaxf/(1+Hm/K_iH);

deno=(1+(A/K_A)^n)*(1+B/K_B+D/K_D)*(1+C/K_C)*(1+E/K_E);  
J_IDH =Vmaxf_prime/(K_A)^n/K_B*((A)^n*B-C*D*E/Keq)/deno;
end