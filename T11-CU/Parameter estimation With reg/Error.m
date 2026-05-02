function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(:,1)/10^6; %Flux(mmol/min)
A = Data(:,2)/10^6;%Ca
C = Data(:,3)/10^3;%Mg
E = Data(:,4)/10^3;%Pi
dPsi=Data(:,5);

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), C(j), E(j),dPsi(j)];
    Ymodel(j,:) = Model(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
