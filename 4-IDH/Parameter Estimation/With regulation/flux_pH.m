function J_CITS=flux_pH(mpar,pH)

Hm = 10^(-pH);                   
%%% Binding constants for 5 parameters
K_H=mpar(1); %proton binding constant(M)
Vmaxf=mpar(2); %max forward velocity

J_CITS=Vmaxf/(1+(Hm/K_H));
end
