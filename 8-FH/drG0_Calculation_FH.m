% Calculation of transformed drG0 of FH biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH_vec = 7;
 
pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The FH reference reaction: 
% FUM(2-) + H2O(0) = MAL(2-)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_FUM = -311.04; dfG0_H2O = -155.66; dfG0_MAL = -682.85;
 
% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (dfG0_MAL) - (dfG0_FUM + dfG0_H2O);
Keq = exp(-drG0/(R*T)); % Keq of reference reaction

disp('Values of drG0 and Keq0 at pH=7:');
disp([drG0, Keq]); 

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
pK_HFUM = 4.09;      % HFUM(-)  = H(+) + FUM(2-)
pK_KFUM = -inf;      % KFUM(-)  = K(+) + FUM(2-)
pK_MgFUM = -inf;     % MgFUM(0) = Mg(2+) + FUM(2-)
pK_CaFUM = 0.60;     % CaFUM(0) = Ca(2+) + FUM(2-)
pK_NaFUM = -inf;     % NaFUM(-) = Na(+) + FUM(2-)

pK_HMAL = 4.71;      % HMAL(-)  = H(+) + MAL(2-)
pK_KMAL = 0.18;      % KMAL(-)  = K(+) + MAL(2-)
pK_MgMAL = 1.71;     % MgMAL(0) = Mg(2+) + MAL(2-)
pK_CaMAL = 1.95;     % CaMAL(0) = Ca(2+) + MAL(2-)
pK_NaMAL = 0.28;     % NaMAL(-) = Na(+) + MAL(2-)

% Free ion concentrations in the solution mixture
K_free = 140e-3; 
Mg_free = 0.8e-3;
Ca_free = 100e-9; 
Na_free = 10e-3;
 
% Account for binding polynomials in thermodynamic parameters calculations
P_FUM = 1 + H_free*10^pK_HFUM + K_free*10^pK_KFUM + Mg_free*10^pK_MgFUM ...
      + Ca_free*10^pK_CaFUM+ Na_free*10^pK_NaFUM;
P_MAL = 1 + H_free*10^pK_HMAL + K_free*10^pK_KMAL + Mg_free*10^pK_MgMAL ...
      + Ca_free*10^pK_CaMAL+ Na_free*10^pK_NaMAL;


Keqp = Keq*(P_MAL/(P_FUM));
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp with pH and ion binding effect:');
disp([drG0p, Keqp]);

