% Calculation of transformed drG0 of SCAS biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH_vec = 7;

pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The SCAS reference reaction: 
% GDP(ADP)(3-) + SCOA(-) + PI(2-) = COA(-) + SUC(2-) + GTP(ATP)(4-) + H(+)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_GDP = -1424.70; dfG0_SCOA = -347.47; dfG0_PI = -1059.49;
dfG0_COA = -7.26; dfG0_SUC = -530.64; dfG0_GTP = -2292.50;
 
% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (dfG0_COA + dfG0_SUC + dfG0_GTP) - ...
    (dfG0_GDP + dfG0_SCOA + dfG0_PI); % kJ/mol (pH=0)
Keq = exp(-drG0/(R*T)); % Keq of reference reaction without pH effect

disp('Values of drG0 and Keq0 at pH =7:');
disp([drG0, Keq]); 

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
pK_HSCOA = 3.99;      % HSCoA(0)  = H(+) + SCoA(-)
pK_KSCOA = -inf;      % KSCoA(0)  = K(+) + SCoA(-)
pK_MgSCOA = -inf;     % MgSCoA(+) = Mg(2+) + SCoA(-)
pK_CaSCOA = -inf;     % CaSCoA(+) = Ca(2+) + SCoA(-)
pK_NaSCOA = -inf;     % NaSCoA(0) = Na(+) + SCoA(-)

pK_HPI = 6.78;      % HPI(-)  = H(+) + PI(2-)
pK_KPI = 0.50;      % KPI(-)  = K(+) + PI(2-)
pK_MgPI = 1.82;     % MgPI(0) = Mg(2+) + PI(2-)
pK_CaPI = 1.74;     % CaPI(0) = Ca(2+) + PI(2-)
pK_NaPI = 0.61;     % NaPI(-) = Na(+) + PI(2-)

pK_HCOA = 8.17;      % HCOA(0)  = H(+) + COA(-)
pK_KCOA = -inf;      % KCOA(0)  = K(+) + COA(-)
pK_MgCOA = -inf;     % MgCOA(+)  = Mg(2+) + COA(-)
pK_CaCOA = -inf;     % CaCOA(+)  = Ca(2+) + COA(-)
pK_NaCOA = -inf;     % NaCOA(0)  = Na(+) + COA(-)

pK_HSUC = 5.27;      % HSUC(-)  = H(+) + SUC(2-)
pK_KSUC = 0.43;      % KSUC(-)  = K(+) + SUC(2-)
pK_MgSUC = 1.35;     % MgSUC(0) = Mg(2+) + SUC(2-)
pK_CaSUC = 1.40;     % CaSUC(0) = Ca(2+) + SUC(2-)
pK_NaSUC = 0.42;     % NaSUC(-) = Na(+) + SUC(2-)

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
P_SCOA = 1 + H_free*10^pK_HSCOA + K_free*10^pK_KSCOA + Mg_free*10^pK_MgSCOA ...
       + Ca_free*10^pK_CaSCOA+ Na_free*10^pK_NaSCOA;
P_PI = 1 + H_free*10^pK_HPI + K_free*10^pK_KPI + Mg_free*10^pK_MgPI ...
       + Ca_free*10^pK_CaPI+ Na_free*10^pK_NaPI; 
P_COA = 1 + H_free*10^pK_HCOA + K_free*10^pK_KCOA + Mg_free*10^pK_MgCOA ...
      + Ca_free*10^pK_CaCOA+ Na_free*10^pK_NaCOA;
P_SUC = 1 + H_free*10^pK_HSUC + K_free*10^pK_KSUC + Mg_free*10^pK_MgSUC ...
      + Ca_free*10^pK_CaSUC+ Na_free*10^pK_NaSUC;
P_ATP = 1 + H_free*10^pK_HATP + K_free*10^pK_KATP + Mg_free*10^pK_MgATP ...
      + Ca_free*10^pK_CaATP+ Na_free*10^pK_NaATP;
P_ADP = 1 + H_free*10^pK_HADP + K_free*10^pK_KADP + Mg_free*10^pK_MgADP ...
      + Ca_free*10^pK_CaADP+ Na_free*10^pK_NaADP;  
  
Keqp = Keq*(P_COA*P_SUC*P_ATP/(P_ADP*P_SCOA*P_PI));
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp with ion binding effect:');
disp([drG0p, Keqp]);

