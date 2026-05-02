function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(:,1)/10^6; %Flux(mmol/min)
A = Data(:,2)/10^3;%
B = Data(:,3)/10^3;%

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j),0];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
