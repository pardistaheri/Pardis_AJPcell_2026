function J_CITS=flux_pH(mpar,pH)

H = 10^(-pH);                   
%%% Binding constants 
K_H=mpar(1); %1st proton binding constant(M) 

J_CITS=1/(1+H./K_H);