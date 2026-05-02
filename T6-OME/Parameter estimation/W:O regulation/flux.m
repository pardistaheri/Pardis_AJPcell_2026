function J_OME=flux(mpar,Conc)

R = 8.314/1000; %gas constant (KJ/K/mol)
T = 310.15; %tempreture in the isolated mito experiment(K)
F  = 0.096484; %Faraday constant(KJ/mol/mV)
RT = R*T;

%%% Binding constants
K_A=mpar(1);  
K_B=mpar(2);
K_C=K_A;
K_D=K_B;
Tmaxf = mpar(3);

%%% metaboloite concentration in the experiment cell
A=Conc(1); %OXA_e(M)
B=Conc(2); %MAL_m(M)
C = Conc(3); %OXA_m(M)
D = Conc(4); %MAL_e(M)

deno = 1+(A/K_A)+(B/K_B)+(C/K_C)+(D/K_D)+(A*B/K_A/K_B)+(C*D/K_C/K_D);
J_OME = Tmaxf * ((A*B/K_A/K_B)-(C*D/K_C/K_D))/deno;
end
