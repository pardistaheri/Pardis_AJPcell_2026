function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(:,1)/10^3; %Flux(mM/min)
A = Data(:,2)/10^3; %C_MALm 
B = Data(:,3)/10^3; %C_NADm 
C = Data(:,4)/10^3; %C_OXAm 
D = Data(:,5)/10^3; %C_NADHm 

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j),C(j),D(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end