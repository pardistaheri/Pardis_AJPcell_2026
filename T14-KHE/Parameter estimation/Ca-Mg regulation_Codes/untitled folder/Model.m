function Tmax_reg=Model(mpar,Conc)

Mg = Conc;

K_Mg=mpar;    

Tmax_reg = 1/(1+(Mg/K_Mg));
end
