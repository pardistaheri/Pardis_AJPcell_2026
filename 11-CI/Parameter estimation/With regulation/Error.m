function [Err] = Error(mpar)

% Load experimental data
Data = load('Data_Ca.txt');
Flux = Data(:,1)/10^6; %Flux(mM/min)
A = 20/10^6; %C_NADH (M)
Ca = Data(:,3)/10^6; %C_Ca (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [Ca(j),A];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
