function J_CITS=flux_pH(mpar,pH)

Hm = 10^(-pH);                   
%%% Binding constants for 5 parameters
K_H=mpar(1); %1st proton binding constant(M)

J_CITS=1/(1+(Hm/K_H));
end
