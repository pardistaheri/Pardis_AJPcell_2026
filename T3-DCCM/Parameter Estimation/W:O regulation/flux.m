function J_DCC=flux(mpar,Conc)

%%% Binding constants
K_A=mpar(1);  
K_B=mpar(2);
K_C=K_A;
K_D=K_B;
Tmaxf = mpar(3);

%%% metaboloite concentration in the experiment cell
A=Conc(1); %Pi_m(M)
B=Conc(2); %MAL_e(M)
C=Conc(3); %Pi_e(M)
D=Conc(4); %MAL_m(M) 

Tmax_Prime = Tmaxf/(1+(C/K_C));
deno = 1+(A/K_A)+(B/K_B)+(C/K_C)+(D/K_D)+(A*B/K_A/K_B)+(C*D/K_C/K_D);
J_DCC = Tmax_Prime * (((A*B/K_A/K_B))-((C*D/K_C/K_D)))/deno;
end
