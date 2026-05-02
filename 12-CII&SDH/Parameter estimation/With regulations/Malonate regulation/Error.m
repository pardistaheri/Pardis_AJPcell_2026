function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(5:10,1); %Flux(mM/min)
A = Data(5:10,2)/10^3; %C_Sucm (M)
MA = Data(5:10,4)/10^6; %C_MAmm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), MA(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end

