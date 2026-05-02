% Reaction 7: Nucleoside diphosphokinase(NDK)- Enzyme
% Nucleoside diphosphokinase flux equation regulated by pH
% GTP + ADP = GDP + ATP 
% A - GTP; B - ADP; ; C - GDP; D - ATP;
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in
% Biomedical Engineering

function J_NDK=NDK(GTPm,ADPm,GDPm,ATPm,pH_m,p)
	
%%% Thermodynamics
dGr= 0;  % kJ/mol Gibbs free energy of the reaction at pH=7  
Keq0=exp(-dGr/(p.R_con*p.Tem));
Keq=Keq0*10^(pH_m-7); % pH correction 

%%% Km parameters 
KA=3.22e-04; % M GTP
KB=1.78e-04; % M ADP
KC=3.22e-04; % M GDP
KD=1.78e-04; % M ATP

%%% Regulation binding constant
KH=1.07e-9; % M H+

%%% Concentration
A=GTPm;
B=ADPm;
C=GDPm;
D=ATPm;
Hm=10^(-pH_m);

%%% flux 

Vmaxf_prime=1/(1+(Hm/KH));

deno=(1+A/KA+C/KC)*(1+B/KB+D/KD);
J_NDK =(Vmaxf_prime/(KA*KB))*(A*B-C*D/Keq)/deno;

