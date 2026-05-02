function J_NHE=Model(mpar,Conc)
% Nam + He -> Nae + Hm

%%% Binding constants
K_A=mpar(1);  
K_B=mpar(2);
K_C = K_A; 
K_D = K_B;
Tmax=mpar(3); %Max forward reaction speed (mmol/min)  
%K_iH =mpar(4);
%nh = 3;


%%% metaboloite concentration in the experiment cell
A=0; %Na_m(M)
B=Conc(1); %H_e(M)
C=Conc(2); %Na_e(M)
D=Conc(3); %H_m(M)

%Tmax_reg = Tmax*(D.^nh/(K_D.^nh+D.^nh));
deno = (1+A/K_A+C/K_C)+(1+B/K_B+D/K_D);
J_NHE = - Tmax*(A*B/K_A/K_B-C*D/K_C/K_D)/deno;
end
