% Calculation of transformed drG0 of ACON biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % temperature (K)
pH_vec = 7;

pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The ACON reference reaction: 
% CIT(3-) = ICIT(3-) 
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_CIT = -966.23; dfG0_ICIT = -959.58; 

% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG = (dfG0_ICIT) - (dfG0_CIT); % kJ/mol (pH=0)
Keq = exp(-drG/(R*T)); % Keq of reference reaction without pH effect

disp('Values of drG0 and Keq0 at pH=7:');
disp([drG, Keq]);  

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
pK_HCIT = 5.67;      % HCIT(2-)  = H(+) + CIT(3-)
pK_KCIT = 0.60;      % KCIT(2-)  = K(+) + CIT(3-)
pK_MgCIT = 3.52;     % MgCIT(-) = Mg(2+) + CIT(3-)
pK_CaCIT = 3.95;     % CaCIT(-) = Ca(2+) + CIT(3-)
pK_NaCIT = 0.75;     % NaCIT(2-) = Na(+) + CIT(3-)

pK_HICIT = 5.76;      % HICIT(2-)  = H(+) + ICIT(3-)
pK_KICIT = -inf;      % KICIT(2-)  = K(+) + ICIT(3-)
pK_MgICIT = 2.62;     % MgICIT(-) = Mg(2+) + ICIT(3-)
pK_CaICIT = 2.54;     % CaICIT(-) = Ca(2+) + ICIT(3-)
pK_NaICIT = -inf;     % NaICIT(2-) = Na(+) + ICIT(3-)

% Free ion concentrations in the solution mixture
K_free = 140e-3; 
Mg_free = 0.8e-3;
Ca_free = 100e-9; 
Na_free = 10e-3;
 
% Account for binding polynomials in thermodynamic parameters calculations
P_CIT = 1 + H_free*10^pK_HCIT + K_free*10^pK_KCIT + Mg_free*10^pK_MgCIT ...
    + Ca_free*10^pK_CaCIT+ Na_free*10^pK_NaCIT;
P_ICIT = 1 + H_free*10^pK_HICIT + K_free*10^pK_KICIT + Mg_free*10^pK_MgICIT ...
    + Ca_free*10^pK_CaICIT+ Na_free*10^pK_NaICIT;

Keqp = Keq*(P_ICIT/(P_CIT));
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp with ion binding effect:');
disp([drG0p, Keqp]);


