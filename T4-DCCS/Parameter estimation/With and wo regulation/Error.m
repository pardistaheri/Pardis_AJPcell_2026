function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(7:16,1); %Flux(mmol/min)
B = Data(7:16,2)/10^3;%
E = Data(7:16,3)/10^3;%

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [B(j),E(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
