function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(:,1); %Flux(%Tmax)
Ca_data = Data(:,2)/10^6; 
Mg_data = Data(:,3)/10^6;

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Ca = Ca_data(j);
    Mg = Mg_data(j);
    Conc = [Ca, Mg];
    Ymodel(j,:) = Model(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
