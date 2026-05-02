function [Err] = Error_F(mpar)

% Load experimental data
Data = load('Data_F.txt');
Flux = Data(:,1)*10^-3; %Flux(mM/min)
A = Data(:,2)/1000; %C_FUMm 
B = Data(:,3)/1000; %C_MALm 

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j)];
    Ymodel(j,:) = Flux_F(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end