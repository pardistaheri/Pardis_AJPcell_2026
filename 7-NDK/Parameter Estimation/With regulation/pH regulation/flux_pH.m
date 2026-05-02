function J_NDK=flux_pH(mpar,pH)

Hm = 10^(-pH);                   
%%% Binding constants for 5 parameters
K_H=mpar(1); %proton binding constant(M)

J_NDK=1/(1+(Hm/K_H));
end
