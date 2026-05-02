function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(:,1)/10^3; %Flux(mM/min)
D = Data(:,2)/10^3; %C_NADHm 
ATP = Data(:,3)/10^3;
ADP = Data(:,4)/10^3;
AMP = Data(:,5)/10^3;

%calculate the model flux based on free concentration
for j=1:length(Flux)
   Conc = [ D(j),ATP(j),ADP(j),AMP(j)];
    Ymodel(j,:) = flux(mpar,Conc); 
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end