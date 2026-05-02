% Calculation of transformed drG0 of CIV biochemical reaction
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
 
% The CIV reference reaction: 
% 2CytCr + 0.5O2 + 2H(+) =  2CytCo + H2O + 2dH
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_CYTC_O = -7.29; dfG0_H2O = -155.66; dfG0_O2 = 16.40; dfG0_CYTC_R = -27.75; 

% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (2*dfG0_CYTC_O + dfG0_H2O)-((2*dfG0_CYTC_R) + (0.5*dfG0_O2)); % kJ/mol 
Keq = exp(-drG0 - 2*F*dPsi/(R*T)); % Keq of reference reaction 

disp('Values of drG0 and Keq at pH = 7:');
disp([drG0, Keq]); 

