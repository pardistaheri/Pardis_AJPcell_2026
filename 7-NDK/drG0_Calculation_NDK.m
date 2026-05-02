% Calculation of transformed drG0 of NDK biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH_vec = 7;
 
pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The NDK reference reaction: 
% GTP(4-) + ADP(3-) = GDP(3-) + ATP(4-)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 0) 
dfG0_ADP = -1903.96; dfG0_ATP = -2771.0;
dfG0_GDP = -1903.96; dfG0_GTP = -2771.0;

% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (dfG0_GDP + dfG0_ATP) - ...
    (dfG0_ADP + dfG0_GTP); % kJ/mol (pH=0)
Keq = exp(-drG0/(R*T)); % Keq of reference reaction without pH effect

disp('Values of drG0 and Keq at pH=0:');
disp([drG0, Keq]); 

Keqp = Keq/H_free; % Keq of reference reaction with pH effect (apparent Keq)
drG0p = -R*T*log(Keqp); % Transformed gibbs free energy of reference reaction

disp('Values of drG0p and Keqp at pH=7:');
disp([drG0p, Keqp]); 

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)

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

pK_HGTP = 6.71;      % HGTP(3-)  = H(+) + GTP(4-)
pK_KGTP = 1.17;      % KGTP(3-)  = K(+) + GTP(4-)
pK_MgGTP = 4.28;     % MgGTP(2-) = Mg(2+) + GTP(4-)
pK_CaGTP = 3.95;     % CaGTP(2-) = Ca(2+) +GTP(4-)
pK_NaGTP = 0.75;     % NaGTP(3-) = Na(+) + GTP(4-)

pK_HGDP = 6.50;      % HGDP(2-)  = H(+) + GDP(3-)
pK_KGDP = 1.00;      % KGDP(2-)  = K(+) + GDP(3-)
pK_MgGDP = 3.30;     % MgGDP(-) = Mg(2+) + GDP(3-)
pK_CaGDP = 2.86;     % CaGDP(-) = Ca(2+) + GDP(3-)
pK_NaGDP = 1.12;     % NaGDP(2-) = Na(+) + GDP(3-)

% Free ion concentrations in the solution mixture
K_free = 140e-3; 
Mg_free = 1e-3;
Ca_free = 500e-9; 
Na_free = 10e-3;
 
% Account for binding polynomials in thermodynamic parameters calculations

P_ATP = 1 + H_free*10^pK_HATP + K_free*10^pK_KATP + Mg_free*10^pK_MgATP ...
      + Ca_free*10^pK_CaATP+ Na_free*10^pK_NaATP;
P_ADP = 1 + H_free*10^pK_HADP + K_free*10^pK_KADP + Mg_free*10^pK_MgADP ...
      + Ca_free*10^pK_CaADP+ Na_free*10^pK_NaADP;  
P_GTP = 1 + H_free*10^pK_HGTP + K_free*10^pK_KGTP + Mg_free*10^pK_MgGTP ...
      + Ca_free*10^pK_CaGTP+ Na_free*10^pK_NaGTP;
P_GDP = 1 + H_free*10^pK_HGDP + K_free*10^pK_KGDP + Mg_free*10^pK_MgGDP ...
      + Ca_free*10^pK_CaGDP+ Na_free*10^pK_NaGDP;    
  
Keqp = Keq*(P_GDP*P_ATP/(P_ADP*P_GTP));
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp:');
disp([drG0p, Keqp]);

