function Vi_V0=flux(mpar,pH)
H = 10^(-pH);
K_H = mpar(1); %Proton binding constant in isolated CII(M)
Vi_V0 = 1/(1+H/K_H);
end
