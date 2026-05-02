% Calculation of transformed drG0 of IDH biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH_vec = 7;
 
pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The IDH reference reaction: 
% ICIT(3-) + NAD(+) + H2O(0) = AKG(2-) + NADH(0) + CO3(2-) + 2H(+)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_ICIT = -959.58; dfG0_NAD = 1059.11; dfG0_H2O = -155.66;
dfG0_AKG = -633.59; dfG0_NADH = 1120.09; dfG0_CO2 = -547.10;
 
% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (dfG0_AKG + dfG0_NADH + dfG0_CO2) - ...
    (dfG0_ICIT + dfG0_NAD + dfG0_H2O); % kJ/mol (pH=0)
Keq = exp(-drG0/(R*T)); % Keq of reference reaction without pH effect

disp('Values of drG0 and Keq0 at pH=7:');
disp([drG0, Keq]); 

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
pK_HICIT = 5.76;      % HICIT(2-)  = H(+) + ICIT(3-)
pK_KICIT = -inf;      % KICIT(2-)  = K(+) + ICIT(3-)
pK_MgICIT = 2.62;     % MgICIT(-) = Mg(2+) + ICIT(3-)
pK_CaICIT = 2.54;     % CaICIT(-) = Ca(2+) + ICIT(3-)
pK_NaICIT = -inf;     % NaICIT(2-) = Na(+) + ICIT(3-)

K_CO2H = 2.71e-3;    % CO2 + H2O = H2CO3 (37 0C, I = 0.25 M)
K_H2CO3 = 2.44e-4;   % H2CO3 = HCO3(-) + H(+) (37 0C, I = 0.25 M)
K_HCO3 = 1.6e-10;    % HCO3(-) = CO3(2-) + H(+) (37 0C, I = 0.25 M)

% Free ion concentrations in the solution mixture
K_free = 140e-3; 
Mg_free = 0.8e-3;
Ca_free = 100e-9; 
Na_free = 10e-3;
 
% Account for binding polynomials in thermodynamic parameters calculations
P_ICIT = 1 + H_free*10^pK_HICIT + K_free*10^pK_KICIT + Mg_free*10^pK_MgICIT ...
         + Ca_free*10^pK_CaICIT+ Na_free*10^pK_NaICIT;
P_NAD = 1; 
P_NADH = 1;
P_AKG = 1;
P_CO2 = 1 + H_free/K_HCO3 + H_free^2/(K_HCO3*K_H2CO3) ...
    + H_free^2/(K_HCO3*K_H2CO3*K_CO2H);

Keqp = Keq*(P_AKG*P_NADH*P_CO2/(P_ICIT*P_NAD));
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp with ion binding effect:');
disp([drG0p, Keqp]);

