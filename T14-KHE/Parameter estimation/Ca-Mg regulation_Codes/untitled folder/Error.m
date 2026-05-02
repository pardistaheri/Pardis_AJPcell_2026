function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(7:13,1); %Flux(%Tmax)
Mg_data = Data(7:13,3)/10^6;

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Mg = Mg_data(j);
    Conc = [Mg];
    Ymodel(j,:) = Model(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
