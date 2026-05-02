% Calculation of transformed drG0 of GOT biochemical reaction
close all; clear all; clc
 
% Standard thermodynamic parameters (fixed)
R = 8.314e-3; % gas constant (kJ/K/mol)
T = 310.15; % tempreture (K)
pH_vec = 7;
 
pH = pH_vec; % pH of reaction mixture
H_free = 10^-pH; % Free H+ concentration
 
% The GOT reference reaction: 
% ASP(2-) + AKG(2-) = GLU(2-) + OXA(2-)
% Gibbs free energy of formation of reference species (kJ/mol) (T = 298.15 K;
% standard reactant concentrations = 1 M, I = 0.17 M, P = 1 atm; pH = 7) 
dfG0_ASP = -452.09; dfG0_AKG = -633.59; dfG0_GLU = -372.15; dfG0_OXA = -715; 

% Gibbs free energy of reference reaction at standard conditions (kJ/mol)
drG0 = (dfG0_GLU + dfG0_OXA) - (dfG0_ASP + dfG0_AKG);
Keq = exp(-drG0/(R*T)); % Keq of reference reaction without pH effect

disp('Values of drG0 and Keq0 at pH =7:');
disp([drG0, Keq]); 

% Binding constants in the dissociation reactions (T = 298 oK and I = 0.17 M)
pK_HASP = 3.65;      % HASP(-)  = H(+) + ASP(2-)
pK_KASP = -inf;      % KASP(-)  = K(+) + ASP(2-)
pK_MgASP = -inf;     % MgASP(0) = Mg(2+) + ASP(2-)
pK_CaASP = -inf;     % CaASP(0) = Ca(2+) + ASP(2-)
pK_NaASP = -inf;     % NaASP(-) = Na(+) + ASP(2-)

pK_HOXA = 3.90;      % HOXA(-)  = H(+) + OXA(2-)
pK_KOXA = -inf;      % KOXA(-)  = K(+) + OXA(2-)
pK_MgOXA = 1.02;     % MgOXA(0) = Mg(2+) + OXA(2-)
pK_CaOXA = 1.60;     % CaOXA(0) = Ca(2+) + OXA(2-)
pK_NaOXA = -inf;     % NaOXA(-) = Na(+) + OXA(2-)

pK_HGLU = 4.06;      % HGLU(0)  = H(+) + GLU(-)
pK_KGLU = -inf;      % KGLU(0)  = K(+) + GLU(-)
pK_MgGLU = 1.82;     % MgGLU(+) = Mg(2+) + GLU(-)
pK_CaGLU = -inf;     % CaGLU(+) = Ca(2+) + GLU(-)
pK_NaGLU = -inf;     % NaGLU(0) = Na(+) + GLU(-)

% Free ion concentrations in the solution mixture
K_free = 140e-3; 
Mg_free = 0.8e-3;
Ca_free = 100e-9; 
Na_free = 10e-3;
 
% Account for binding polynomials in thermodynamic parameters calculations
P_ASP = 1 + H_free*10^pK_HASP + K_free*10^pK_KASP + Mg_free*10^pK_MgASP ...
    + Ca_free*10^pK_CaASP+ Na_free*10^pK_NaASP;
P_OXA = 1 + H_free*10^pK_HOXA + K_free*10^pK_KOXA + Mg_free*10^pK_MgOXA ...
      + Ca_free*10^pK_CaOXA+ Na_free*10^pK_NaOXA;
P_GLU = 1 + H_free*10^pK_HGLU + K_free*10^pK_KGLU + Mg_free*10^pK_MgGLU ...
    + Ca_free*10^pK_CaGLU+ Na_free*10^pK_NaGLU;
P_AKG = 1;

Keqp = Keq*(P_OXA*P_GLU/(P_ASP*P_AKG));
drG0p = -R*T*log(Keqp);

disp('Values of drG0p and Keqp with pH and ion binding effect:');
disp([drG0p, Keqp]);
