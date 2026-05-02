function J_DCC=flux(mpar,Conc)

%%% Binding constants

K_A=3.3280506e-03; %Pim Binding constant (M)(from DCCM estimation)
K_B = 4.9299098e-04; %SUCe Binding constant (M) (from DCCM SUC regulation estimation)
K_C=K_A; %Pie Binding constant (M)
K_D=K_B;
K_E = 2.0911990e-04; %Malonate Binding constant (M)
Tmaxf = mpar(1);

%%% metaboloite concentration in the experiment cell
A=15/1000; %Pi_m(M)
B=Conc(1); %Suc_e(M)
C=0; %Pi_e(M)
D=0; %Suc_m(M) 
E=Conc(2);

K_B_Prime = K_B*(1+(E/K_E));
deno = 1+(A/K_A)+(B/K_B_Prime)+(C/K_C)+(D/K_D)+(A*B/K_A/K_B_Prime)+(C*D/K_C/K_D);
J_DCC = Tmaxf * (((A*B/K_A/K_B_Prime))-((C*D/K_C/K_D)))/deno;
end
