% Calculation of transformed drG0 of CITS biochemical reaction
close all; clear all; clc

% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH_vec = 7; 

pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The CITS reference reaction: 
% ACOA(0) + OXA(2-) + H2O(0) = COA(-) + CIT(3-) + 2H(+)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_ACOA = -66.33; dfG0_OXA = -715; dfG0_H2O = -155.66;
dfG0_COA = -7.26; dfG0_CIT = -966.23;
 
% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG = (dfG0_COA + dfG0_CIT) - ...
    (dfG0_ACOA + dfG0_OXA + dfG0_H2O);% kJ/mol (pH=0)
Keq = exp(-drG/(R*T)); % Keq of reference reaction without pH effect

disp('Values of drG0 and Keq0 wat pH=7:');
disp([drG, Keq]);

% Binding constants in the dissociation reactions 
pK_HOXA = 3.90;      % HOXA (-)  = H(+) + OXA (2-)
pK_KOXA = -inf;      % KOXA (-)  = K(+) + OXA (2-)
pK_MgOXA = 1.02;     % MgOXA (0) = Mg(2+) + OXA (2-)
pK_CaOXA = 1.60;     % CaOXA (0) = Ca(2+) + OXA (2-)
pK_NaOXA = -inf;     % NaOXA (-) = Na(+) + OXA (2-)
 
pK_HCOA = 8.17;      % HCOA(0)  = H(+) + COA(-)
pK_KCOA = -inf;      % KCOA(0)  = K(+) + COA(-)
pK_MgCOA = -inf;     % MgCOA(+)  = Mg(2+) + COA(-)
pK_CaCOA = -inf;     % CaCOA(+)  = Ca(2+) + COA(-)
pK_NaCOA = -inf;     % NaCOA(0)  = Na(+) + COA(-)

pK_HCIT = 5.67;      % HCIT (2-)  = H(+) + CIT (3-)
pK_KCIT = 0.6;       % KCIT (2-)  = K(+) + CIT (3-)
pK_MgCIT = 3.52;     % MgCIT (-)  = Mg(2+) + CIT (3-)
pK_CaCIT = 3.95;     % CaCIT (-)  = Ca(2+) + CIT (3-)
pK_NaCIT = 0.75;     % NaCIT (2-)  = Na(+) + CIT (3-)
 
% Free ion concentrations in the solution mixture
K_free = 140e-3; 
Mg_free = 0.8e-3;
Ca_free = 100e-9; 
Na_free = 10e-3;
 
% Account for binding polynomials in thermodynamic parameters calculations
P_OXA = 1 + H_free*10^pK_HOXA + K_free*10^pK_KOXA + Mg_free*10^pK_MgOXA...
    + Ca_free*10^pK_CaOXA+ Na_free*10^pK_NaOXA;
P_COA = 1 + H_free*10^pK_HCOA + K_free*10^pK_KCOA + Mg_free*10^pK_MgCOA ...
    + Ca_free*10^pK_CaCOA+ Na_free*10^pK_NaCOA;
P_CIT = 1 + H_free*10^pK_HCIT + K_free*10^pK_KCIT + Mg_free*10^pK_MgCIT...
    + Ca_free*10^pK_CaCIT+ Na_free*10^pK_NaCIT;
P_ACOA = 1;

Keqp = Keq*(P_COA*P_CIT/(P_ACOA*P_OXA)); %/(H_free*10^pK_HCOA);
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp with ion binding effect:');
disp([drG0p, Keqp]);


