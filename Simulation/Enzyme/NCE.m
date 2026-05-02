% Transport 15- NCE (Full Version) 
% Calcium-Sodium antiporter with complete proton binding terms
% Cam + 3Nae -> Cae + 3Nam
% Based on Equations B45-B49 from the attached image
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function T_NCE = NCE(Cam,Nae,Cae,Nam,dPsi,Hm,He,p)

RT = p.R_con * p.Tem;

%%% Kinetic parameters from Table T13
Tmax_NCE = 7.84e-3;      % Maximum forward transport rate (mol/s.mg protein)
KA = 1.1e-7;             % Dissociation constant of Ca2+ for antiporter (M)
KC = 1.1e-7;             % Dissociation constant of Ca2+ for antiporter (M)
KB = 2.2e-3;             % Dissociation constant of Na+ for antiporter (M)
KD = 2.2e-3;             % Dissociation constant of Na+ for antiporter (M)
KH1 = 6.4e-8;            % 1st proton binding constant (M)
KH2 = 1.4e-7;            % 2nd proton binding constant (M)

%%% Concentrations
A = Cam;    % Ca2+ mitochondrial matrix (M)
B = Nae;    % Na+ extracellular (M)
C = Cae;    % Ca2+ extracellular (M)
D = Nam;    % Na+ mitochondrial matrix (M)

%%% Modified dissociation constants
K_A_prime = KA * ( (Hm/KH1) + (Hm^2/KH1^2) + (Hm^3/KH1^3) + 1 + (KH2/Hm) + (KH2^2/Hm^2) + (KH2^3/Hm^3) );
K_C_prime = KC * ( (He/KH1) + (He^2/KH1^2) + (He^3/KH1^3) + 1 + (KH2/He) + (KH2^2/He^2) + (KH2^3/He^3) );


D_Hm = K_A_prime * ((Hm / KH2)^3); 
D_He = K_C_prime * ((He / KH2)^3); 


% Mitochondrial side denominator (Hm term)
for n = 1:3
    denom_Hm = 1 + (A/K_A_prime)+ (B^n / KB^n) + (A * B^n) / (K_A_prime * KB^n);
end

% Extracellular side denominator (He term)
for n = 1:3
    denom_He = 1 + (C/K_C_prime); + (D^n / KD^n) + (C * D^n) / (K_C_prime * KD^n);
end

denominator = D_Hm * denom_Hm + D_He * denom_He;

%%% Numerator
term1 = Tmax_NCE * exp(0.5 * p.F_con * dPsi / RT) * D_Hm * (A * B^3) / (K_A_prime * KB^3);
term2 = exp(-0.5 * p.F_con * dPsi / RT) * D_He * (C * D^3) / (K_C_prime * KD^3);

numerator = term1 - term2;

%%% Final transport rate
T_NCE = numerator / denominator;
