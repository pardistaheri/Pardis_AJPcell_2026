function J_KHE=Model(mpar,Conc)
%Km + He -> Ke + Hm

%%% Binding constants

K_A = mpar(1);
K_B = mpar(2);
K_C = K_A;  
K_D = K_B;
Tmax=mpar(3); %Max forward reaction speed (mmol/min)  

%%% metaboloite concentration in the experiment cell
A=Conc(1); %K_m(M)
B=Conc(2); %H_e(M)
C=Conc(3); %K_e(M)
D=Conc(4); %H_m(M)

deno = 1+A/K_A+C/K_C+B/K_B+D/K_D+(A*B/(K_A*K_B))+((C*D)/(K_C*K_D));
J_KHE = -Tmax*((A*B/(K_A*K_B))-((C*D)/(K_C*K_D)))/deno;
end
