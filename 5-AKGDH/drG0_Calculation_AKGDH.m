% Calculation of transformed drG0 of AKGDH biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH_vec = 7;
 
pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The AKGDH reference reaction: 
% AKG(2-) + NAD(+) + COA(-) + H2O = SCOA(-) + NADH(0) + CO3(2-) + H(+)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_AKG = -633.59; dfG0_NAD = 1059.11; dfG0_COA = -7.26; dfG0_H2O = -155.66;
 dfG0_SCOA = -347.47; dfG0_NADH = 1120.09; dfG0_CO2 = -547.10;
 
% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (dfG0_SCOA + dfG0_NADH + dfG0_CO2) - ...
    (dfG0_AKG + dfG0_NAD + dfG0_COA + dfG0_H2O); % kJ/mol (pH=0)
Keq = exp(-drG0/(R*T)); % Keq of reference reaction without pH effect

disp('Values of drG0 and Keq0 at pH=7:');
disp([drG0, Keq]); 

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
pK_HCOA = 8.17;      % HCOA(0)  = H(+) + COA(-)
pK_KCOA = -inf;      % KCOA(0)  = K(+) + COA(-)
pK_MgCOA = -inf;     % MgCOA(+)  = Mg(2+) + COA(-)
pK_CaCOA = -inf;     % CaCOA(+)  = Ca(2+) + COA(-)
pK_NaCOA = -inf;     % NaCOA(0)  = Na(+) + COA(-)

pK_HSCOA = 3.99;      % HSCoA(0)  = H(+) + SCoA(-)
pK_KSCOA = -inf;      % KSCoA(0)  = K(+) + SCoA(-)
pK_MgSCOA = -inf;     % MgSCoA(+) = Mg(2+) + SCoA(-)
pK_CaSCOA = -inf;     % CaSCoA(+) = Ca(2+) + SCoA(-)
pK_NaSCOA = -inf;     % NaSCoA(0) = Na(+) + SCoA(-)

K_CO2H = 2.71e-3;    % CO2 + H2O = H2CO3 (37 0C, I = 0.25 M)
K_H2CO3 = 2.44e-4;   % H2CO3 = HCO3(-) + H(+) (37 0C, I = 0.25 M)
K_HCO3 = 1.6e-10;    % HCO3(-) = CO3(2-) + H(+) (37 0C, I = 0.25 M)


% Free ion concentrations in the solution mixture
K_free = 140e-3; 
Mg_free = 0.8e-3;
Ca_free = 100e-9; 
Na_free = 10e-3;
 
% Account for binding polynomials in thermodynamic parameters calculations
P_COA = 1 + H_free*10^pK_HCOA + K_free*10^pK_KCOA + Mg_free*10^pK_MgCOA ...
      + Ca_free*10^pK_CaCOA+ Na_free*10^pK_NaCOA;
P_SCOA = 1 + H_free*10^pK_HSCOA + K_free*10^pK_KSCOA + Mg_free*10^pK_MgSCOA ...
       + Ca_free*10^pK_CaSCOA+ Na_free*10^pK_NaSCOA;
P_NADH = 1;
P_AKG = 1;
P_NAD = 1;
P_CO2 = 1 + H_free/K_HCO3 + H_free^2/(K_HCO3*K_H2CO3) ...
      + H_free^2/(K_HCO3*K_H2CO3*K_CO2H);
 
Keqp = Keq*(P_SCOA*P_NADH*P_CO2/(P_AKG*P_NAD*P_COA));
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp with ion binding effect:');
disp([drG0p, Keqp]);

