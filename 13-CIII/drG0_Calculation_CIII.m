% Calculation of transformed drG0 of CIII biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
dPsi = -140; % (mV)
%dSi = -0.14; % (V)
F = 0.096484; % (kJ/mol/mV)
pH= 7;
 
H_free = 10^-pH; % Free H+ concentration
 
% The CIII reference reaction: 
% UQH2(0) + 2CYTC_O(3+) + 2H(+) = UQ(0) + 2CYTC_R(2+) + 4dH
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_UQH2 = 3660.55; dfG0_CYTC_O = -7.29; dfG0_UQ = 3668.94; dfG0_CYTC_R = -27.75; 

% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (2*dfG0_CYTC_R + dfG0_UQ) - (dfG0_UQH2 + 2*dfG0_CYTC_O); % kJ/mol
Keq0 = exp(-drG0 - 2*F*dPsi/(R*T)); % Keq of reference reaction 

disp('Values of drG0 and Keq:');
disp([drG0, Keq0]); 


