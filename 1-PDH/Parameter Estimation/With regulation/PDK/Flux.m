%flux equation
function PDK=Flux(mpar,C)
A = C(:,1);
ATP = C(:,2);
ADP = C(:,3);
pH = C(:,4);
H=10^(-pH);

% Regulatory Constants for kinase
K_H=mpar(1);
K_ATP=mpar(2);
K_ADP=mpar(3);
K_iPYR=mpar(4);


% kinase activity (normalized)
PDK = (ATP./K_ATP)./((H./K_H+K_H./H).*(1+(A./K_iPYR)).*(1 +(ATP./K_ATP)+(ADP./K_ADP)));
end