%flux equation
function PDP=Flux(Cation, mpar)
H = Cation(1);
Mg = Cation(2);
Ca = Cation(3);

% Regulatory Constants for phosphatase
KMg=mpar(1);
KH=mpar(2);
KCa=mpar(3);
bCa=mpar(4);


% phosphatase activity (normalized)
PDP = ((Mg./KMg) + (bCa.*Ca./KCa))./((H./KH+KH./H).*(1 + (Mg./KMg) + (Ca./KCa)));

end
 