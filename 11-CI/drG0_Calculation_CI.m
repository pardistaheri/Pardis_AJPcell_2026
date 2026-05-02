% Calculation of transformed drG0 of CI biochemical reaction
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
 
% The CI reference reaction: 
% NADH(0) + UQ(0) + H(+) = NAD(+) + UQH2(0) + 4dH
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_NADH = 1120.09; dfG0_UQ = 3668.94; dfG0_NAD = 1059.11; dfG0_UQH2 = 3660.55; 

% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (dfG0_NAD + dfG0_UQH2) - (dfG0_NADH + dfG0_UQ) % kJ/mol 
Keq = exp(-drG0 - 4*F*dPsi/(R*T)); % Keq of reference reaction 

disp('Values of drG0 and Keq0 at pH=7:');
disp([drG0, Keq]); 
