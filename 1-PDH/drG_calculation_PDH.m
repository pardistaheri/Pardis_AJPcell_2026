% Calculation of transformed drG0 of PDH biochemical reaction
close all; clear all; clc
format short e

% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH = 7;

H_free = 10^-pH; % Free H+ concentration

%%% The PDH reference reaction: 
% Pyr(-) + CoA(-) + NAD(-) + H2O(0) = ACoA(0) + CO2(2-) + NADH(2-) + H(+)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_Pyr = -350.78; dfG0_CoA = -7.26; dfG0_NAD = 1059.11; dfG0_H2O = -155.66;
dfG0_ACoA = -66.33; dfG0_CO2 = -547.10; dfG0_NADH = 1120.09; 

% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG = (dfG0_ACoA + dfG0_CO2 + dfG0_NADH) - ...
    (dfG0_Pyr + dfG0_CoA + dfG0_NAD + dfG0_H2O); % 19.59 kJ/mol (pH=0)
Keq = exp(-drG/(R*T)); % Keq of reference reaction pH=7

disp('Values of drG and Keq at pH=7:');
disp([drG, Keq]); 

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
pK_HPyr = 2.26;      % HPyr(0)  = H(+) + Pyr(-)
pK_KPyr = -inf;      % KPyr(0)  = K(+) + Pyr(-)
pK_MgPyr = 1.1;      % MgPyr(+) = Mg(2+) + Pyr(-)
pK_CaPyr = 0.8;      % CaPyr(+) = Ca(2+) + Pyr(-)
pK_NaPyr = -inf;     % NaPyr(0) = Na(+) + Pyr(-)

pK_HCoA = 8.17;      % HCoA(0)  = H(+) + CoA(-)
pK_KCoA = -inf;      % KCoA(0)  = K(+) + CoA(-)
pK_MgCoA = -inf;     % MgCoA(+) = Mg(2+) + CoA(-)
pK_CaCoA = -inf;     % CaCoA(+) = Ca(2+) + CoA(-)
pK_NaCoA = -inf;     % NaCoA(0) = Na(+) + CoA(-)

K_CO2H = 2.71e-3;    % CO2 + H2O = H2CO3 (37 0C, I = 0.25 M)
K_H2CO3 = 2.44e-4;   % H2CO3 = HCO3(-) + H(+) (37 0C, I = 0.25 M)
K_HCO3 = 1.6e-10;    % HCO3(-) = CO3(2-) + H(+) (37 0C, I = 0.25 M)

P_CO2 = 50;                 % Partial pressure of CO2 in cells (mmHg)
S_CO2 = 3.3e-5;             % Solubility of CO2 in cells water (M/mmHg)
C_CO2 = S_CO2*P_CO2;        % Concentration of dissolved CO2 in cells (M) 
C_H2CO3 = K_CO2H*C_CO2;     % Concentration of H2CO3 in cells (M)
C_HCO3 = K_H2CO3*C_H2CO3/H_free; % Concentration of HCO3 in cells (M)
C_CO3 = K_HCO3*C_HCO3/H_free;    % Concentration of CO3 in cells (M)

% disp('Values of C_CO2, C_H2CO3, C_HCO3, and C_CO3:')
% disp([C_CO2, C_H2CO3, C_HCO3, C_CO3]);

% Free ion concentrations in the solution mixture (typical in mitochondria)
K_free = 140e-3; 
Mg_free = 0.8e-3;
Ca_free = 100e-9; 
Na_free = 10e-3;

% Account for binding polynomials in thermodynamic parameters calculations
P_Pyr = 1 + H_free*10^pK_HPyr + K_free*10^pK_KPyr + Mg_free*10^pK_MgPyr ...
    + Ca_free*10^pK_CaPyr+ Na_free*10^pK_NaPyr;
P_COA = 1 + H_free*10^pK_HCoA + K_free*10^pK_KCoA + Mg_free*10^pK_MgCoA ...
    + Ca_free*10^pK_CaCoA+ Na_free*10^pK_NaCoA;
P_CO2 = 1 + H_free/K_HCO3 + H_free^2/(K_HCO3*K_H2CO3) ...
    + H_free^2/(K_HCO3*K_H2CO3*K_CO2H);
P_NAD = 1;
P_ACOA = 1;
P_NADH = 1;

Sum1 = C_CO2 + C_HCO3;
Sum2 = C_CO2 + C_H2CO3 + C_HCO3 + C_CO3;
Sum3 = C_CO3*P_CO2;
Ratio1 = C_HCO3/C_CO3; 
Ratio2 = C_H2CO3/C_CO3; 
Ratio3 = C_CO2/C_CO3;

Keqp = Keq*(P_ACOA*P_NAD*P_CO2/(P_Pyr*P_COA*P_NAD)); %*H_free*10^pK_HCoA;
drG0p = -R*T*log(Keqp);
disp('Values of drG0p and Keqp with ion binding effect:');
disp([drG0p, Keqp]);

disp('Values of CO2 concnetration:')
disp([Sum2]);
