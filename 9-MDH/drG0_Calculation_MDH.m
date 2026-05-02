% Calculation of transformed drG0 of MDH biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH_vec = 7;
 
pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The MDH reference reaction: 
% MAL(2-) + NAD(+) = OXA(2-) + NADH(0) + H(+)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_MAL = -682.85; dfG0_NAD = 1059.11; dfG0_OXA = -715; dfG0_NADH = 1120.09;

% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (dfG0_OXA + dfG0_NADH) - (dfG0_MAL + dfG0_NAD; 
Keq = exp(-drG0/(R*T)); % Keq of reference reaction without pH effect

disp('Values of drG0 and Keq at pH=7:');
disp([drG0, Keq]); 

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
pK_HMAL = 4.71;      % HMAL(-)  = H(+) + MAL(2-)
pK_KMAL = 0.18;      % KMAL(-)  = K(+) + MAL(2-)
pK_MgMAL = 1.71;     % MgMAL(0) = Mg(2+) + MAL(2-)
pK_CaMAL = 1.95;     % CaMAL(0) = Ca(2+) + MAL(2-)
pK_NaMAL = 0.28;     % NaMAL(-) = Na(+) + MAL(2-)

pK_HOXA = 3.90;      % HOXA(-)  = H(+) + OXA(2-)
pK_KOXA = -inf;      % KOXA(-)  = K(+) + OXA(2-)
pK_MgOXA = 1.02;     % MgOXA(0) = Mg(2+) + OXA(2-)
pK_CaOXA = 1.60;     % CaOXA(0) = Ca(2+) + OXA(2-)
pK_NaOXA = -inf;     % NaOXA(-) = Na(+) + OXA(2-)

% Free ion concentrations in the solution mixture
K_free = 140e-3; 
Mg_free = 0.8e-3;
Ca_free = 100e-9; 
Na_free = 10e-3;
 
% Account for binding polynomials in thermodynamic parameters calculations

P_MAL = 1 + H_free*10^pK_HMAL + K_free*10^pK_KMAL + Mg_free*10^pK_MgMAL ...
      + Ca_free*10^pK_CaMAL+ Na_free*10^pK_NaMAL;
P_OXA = 1 + H_free*10^pK_HOXA + K_free*10^pK_KOXA + Mg_free*10^pK_MgOXA ...
      + Ca_free*10^pK_CaOXA+ Na_free*10^pK_NaOXA;
P_NAD = 1; 
P_NADH = 1;

Keqp = Keq*(P_OXA*P_NADH/(P_MAL*P_NAD));
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp with pH and ion binding effect:');
disp([drG0p, Keqp]);

