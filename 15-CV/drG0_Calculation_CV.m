% Calculation of transformed drG0 of CV biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
dPsi = -140; % (mV)
%dSi = -0.14; % (V)
F = 0.096484; % (kJ/mol/mV)
pH_vec = 7;
 
pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The CV reference reaction: 
% ADP + Pi + H2O =  ATP
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_ATP = -2292.50; dfG0_H2O = -155.66; dfG0_ADP = -1424.70; dfG0_Pi = -1059.49; 
%dfG0_ATP = -2771.00; dfG0_H2O = -235.74; dfG0_ADP = -1903.96; dfG0_Pi = -1098.27; 
% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (dfG0_ATP+dfG0_H2O)-(dfG0_Pi + dfG0_ADP); % kJ/mol 
Keq = exp(-drG0 - 2*F*dPsi/(R*T)); % Keq of reference reaction

disp('Values of drG0 and Keq at pH=7:');
disp([drG0, Keq]); 

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
pK_HPI = 6.78;      % HPI(-)  = H(+) + PI(2-)
pK_KPI = 0.50;      % KPI(-)  = K(+) + PI(2-)
pK_MgPI = 1.82;     % MgPI(0) = Mg(2+) + PI(2-)
pK_CaPI = 1.74;     % CaPI(0) = Ca(2+) + PI(2-)
pK_NaPI = 0.61;     % NaPI(-) = Na(+) + PI(2-)

pK_HATP = 6.71;      % HATP(3-)  = H(+) + ATP(4-)
pK_KATP = 1.17;      % KATP(3-)  = K(+) + ATP(4-)
pK_MgATP = 4.28;     % MgATP(2-) = Mg(2+) + ATP(4-)
pK_CaATP = 3.95;     % CaATP(2-) = Ca(2+) + ATP(4-)
pK_NaATP = 0.75;     % NaATP(3-) = Na(+) + ATP(4-)

pK_HADP = 6.50;      % HADP(2-)  = H(+) + ADP(3-)
pK_KADP = 1.00;      % KADP(2-)  = K(+) + ADP(3-)
pK_MgADP = 3.30;     % MgADP(-) = Mg(2+) + ADP(3-)
pK_CaADP = 2.86;     % CaADP(-) = Ca(2+) + ADP(3-)
pK_NaADP = 1.12;     % NaADP(2-) = Na(+) + ADP(3-)

% Free ion concentrations in the solution mixture
K_free = 140e-3; 
Mg_free = 0.8e-3;
Ca_free = 100e-9; 
Na_free = 10e-3;
 
% Account for binding polynomials in thermodynamic parameters calculations
P_PI = 1 + H_free*10^pK_HPI + K_free*10^pK_KPI + Mg_free*10^pK_MgPI ...
       + Ca_free*10^pK_CaPI+ Na_free*10^pK_NaPI; 
P_ATP = 1 + H_free*10^pK_HATP + K_free*10^pK_KATP + Mg_free*10^pK_MgATP ...
      + Ca_free*10^pK_CaATP+ Na_free*10^pK_NaATP;
P_ADP = 1 + H_free*10^pK_HADP + K_free*10^pK_KADP + Mg_free*10^pK_MgADP ...
      + Ca_free*10^pK_CaADP+ Na_free*10^pK_NaADP;  
  
Keqp = Keq*(P_ATP/(P_ADP*P_PI));
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp with ion binding effect:');
disp([drG0p, Keqp]);