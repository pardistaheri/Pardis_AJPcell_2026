function Tmax_reg=Model(mpar,pH)

Hm = 10^(-pH);

K_Hi=mpar(1); %1st proton binding constant(M)
Tmax=mpar(2); %Max forward reaction speed(mmol/min)   
nh =mpar(3);

Tmax_reg = Tmax*(Hm.^nh/(K_Hi.^nh+Hm.^nh));
end
