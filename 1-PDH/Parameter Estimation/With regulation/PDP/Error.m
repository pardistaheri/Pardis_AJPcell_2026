%error function for fmincon
function Err = Error(mpar)
Data = load('Data.txt');
Flux_data = Data(:,1);
pH = Data(:,2);
Mg = Data(:,3)/1000;          
Ca = Data(:,4)/1000000;      

for i=1:1:length(Flux_data)
    %%% Experimental conditions (total cations)
    H_tot = 10^-pH(i);       % pH = 7.5 for Tsai et al. 1973
    Mg_tot = Mg(i);          % Mg_tot = 1 mM for Tsai et al. 1973 (M)
    Ca_tot = Ca(i);          % Ca_tot = 2.6 mM for Tsai et al. 1973
    
    Cation_frac = [H_tot,Mg_tot,Ca_tot];
    Ymodel(i,:) = Flux(Cation_frac, mpar);  
    %calculate the model flux based on free concentration
    
end
Err = sum(((Ymodel-Flux_data)./(Flux_data)).^2);
end