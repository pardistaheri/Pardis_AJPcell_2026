% Calculation of transformed drG0 of CII biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH_vec = 7;

pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The CII/SDH reference reaction: 
% SUC + UQ(0) = FUM(2+) + UQH2(0)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
%dfG0_SUC = -530.64; dfG0_UQ = 3668.94; dfG0_FUM = -311.04; dfG0_UQH2 = 3660.55; 
dfG0_SUC = -690.44; dfG0_UQ = 65.17; dfG0_FUM = -603.32; dfG0_UQH2 = -23.30; 
% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG = (dfG0_FUM + dfG0_UQH2) - (dfG0_SUC + dfG0_UQ) % kJ/mol 
Keq = exp(-drG/(R*T)) % Keq of reference reaction

% disp('Values of drG0 and Keq at pH=7:');
% disp([drG, Keq]); 

Keqp = Keq/H_free % Keq of reference reaction with pH effect (apparent Keq)
drG0p = -R*T*log(Keqp) % Transformed gibbs free energy of reference reaction

% disp('Values of drG0p and Keqp:');
% disp([drG0p, Keqp]); 

% % Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
% pK_HSUC = 5.27;      % HSUC(-)  = H(+) + SUC(2-)
% pK_KSUC = 0.43;      % KSUC(-)  = K(+) + SUC(2-)
% pK_MgSUC = 1.35;     % MgSUC(0) = Mg(2+) + SUC(2-)
% pK_CaSUC = 1.40;     % CaSUC(0) = Ca(2+) + SUC(2-)
% pK_NaSUC = 0.42;     % NaSUC(-) = Na(+) + SUC(2-)
% 
% pK_HFUM = 4.09;      % HATP(3-)  = H(+) + ATP(4-)
% pK_KFUM = -inf;      % KATP(3-)  = K(+) + ATP(4-)
% pK_MgFUM = -inf;     % MgATP(2-) = Mg(2+) + ATP(4-)
% pK_CaFUM = 0.6;     % CaATP(2-) = Ca(2+) + ATP(4-)
% pK_NaFUM = -inf;     % NaATP(3-) = Na(+) + ATP(4-)
% 
% % Free ion concentrations in the solution mixture
% K_free = 140e-3; 
% Mg_free = 0.8e-3;
% Ca_free = 100e-9; 
% Na_free = 10e-3;
%  
% % Account for binding polynomials in thermodynamic parameters calculations
% P_SUC = 1 + H_free*10^pK_HSUC + K_free*10^pK_KSUC + Mg_free*10^pK_MgSUC ...
%       + Ca_free*10^pK_CaSUC+ Na_free*10^pK_NaSUC;
% P_FUM = 1 + H_free*10^pK_HFUM + K_free*10^pK_KFUM + Mg_free*10^pK_MgFUM ...
%       + Ca_free*10^pK_CaFUM+ Na_free*10^pK_NaFUM;
% 
%   
% Keqp = Keq*(P_FUM/(P_SUC));
% drG0p = -R*T*log(Keqp);
% 
% disp('Values of drG0p and Keqp with ion binding effect:');
% disp([drG0p, Keqp]);


