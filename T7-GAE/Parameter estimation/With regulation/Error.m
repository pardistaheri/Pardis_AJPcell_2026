function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(:,1)/10^6; %Flux(mmol/min)
A = Data(:,2)/10^3;%ASP
B = 10/10^3;%Glu
Ca = Data(:,3)/10^6;
%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B, Ca(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
