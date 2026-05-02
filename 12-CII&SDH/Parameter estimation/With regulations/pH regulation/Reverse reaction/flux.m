function Vi_V0=Model(mpar,pH)
H = 10^(-pH);
K_H1 = mpar(1); %Proton binding constant in isolated CII(M)
K_H2 = mpar(2); %Proton binding constant in isolated CII(M)
Vi_V0 = 1/(K_H1/H+1+H/K_H2);
end
