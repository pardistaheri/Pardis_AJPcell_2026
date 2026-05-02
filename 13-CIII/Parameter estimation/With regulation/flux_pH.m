function Vi_V0=flux_pH(mpar,pH)
H = 10^(-pH);
K_H1 = mpar(1); %Proton binding constant in isolated CIII(M)
K_H2 = mpar(2); %Proton binding constant in isolated CIII(M)
Vi_V0 = 1/(K_H1/H+1+H/K_H2);
end
