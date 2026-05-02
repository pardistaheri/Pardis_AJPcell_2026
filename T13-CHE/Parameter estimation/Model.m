function J_CHE=Model(mpar,Conc)
% Cam + He -> Cae + Hm

%%% Binding constants
K_A=mpar(1);  
K_B = 10^-7;
K_C = mpar(2);
K_D = K_B;
Tmax=mpar(3); %Max forward reaction speed (mmol/min)  

%%% metaboloite concentration in the experiment cell
A=Conc(1); %Ca_m(M)
B=10^-7.15; %H_e(M)
C=Conc(2); %Ca_e(M)
D=10^-7.35; %H_m(M)
% 
% J_CHE = Tmax * ((B^2*A - D^2*C) / ...
%         (Km*(B^2+D^2) + B^2*A + D^2*C));
    
%J_CHE = Tmax *(A*B^2 - C*D^2) ./ (K_A*(B^2+D^2) + A*B^2 + C*D^2);
deno = (1+A/K_A+C/K_C)+(1+(B/K_B)^2+(D/K_D)^2);
J_CHE = Tmax*((A*(B^2)/(K_A*(K_B^2)))-((C*(D^2))/(K_C*(K_D^2))))/deno;
end
