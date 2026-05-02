function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(:,1)/10^3; %Flux(mM/min)
A = Data(:,2)/10^3; %C_Sucm (M)
OAA = Data(:,4)/10^6; %C_OAAmm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), OAA(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
%Err = sum(((Cdata-Cmodel)./Cdata).^2);
