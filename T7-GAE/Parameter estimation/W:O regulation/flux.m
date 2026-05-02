function J_GAE=flux(mpar,Conc)

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
A=Conc(1); %ASP_m(M)
B=Conc(2); %GLU_e(M)
C = 0; %ASP_e(M)
D = 0; %GLU_m(M)
dPsi = Conc(3); 
beta = 0.7;

deno = 1+(A/K_A)+(B/K_B)+(C/K_C)+(D/K_D)+(A*B/K_A/K_B)+(C*D/K_C/K_D);
J_GAE = Tmaxf * ((exp(beta*F*dPsi/RT)*(A*B/K_A/K_B))-((exp(-(1-beta)*F*dPsi/RT)*(C*D/K_C/K_D))))/deno;
end
