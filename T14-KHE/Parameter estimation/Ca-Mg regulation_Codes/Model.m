function Tmax_reg=Model(mpar,Conc)

Ca = Conc(1);
Mg = Conc(2);

K_Ca=mpar(1); 
K_Mg=mpar(2);    

Tmax_reg = 1/(1+(Ca/K_Ca)+(Mg/K_Mg));
end
