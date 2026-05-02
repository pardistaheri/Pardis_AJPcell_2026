function J_CU=Model(mpar,Conc)
% Cae -> Cam

%%% Binding constants
K_A=mpar(1);  
K_B=mpar(1);
Tmaxf = mpar(2);

%%% metaboloite concentration in the experiment cell
A=Conc(1); %Ca_e(M)
B=0; %Ca_m(M)

deno = 1+(A^2/K_A^2)+(B^2/K_B^2);
J_CU = Tmaxf * ((A^2/K_A^2)-(B^2/K_B^2))/deno;

end
