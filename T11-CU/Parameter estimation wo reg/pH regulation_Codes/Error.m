function [Err] = Error(mpar)

% Load experimental data
Data = load('Data_pH.txt');
Flux = Data(:,1)/10^6; %Flux(%Tmax)
pH_data = Data(:,2); 

%calculate the model flux based on free concentration
for j=1:length(Flux)
    pH = pH_data(j);
    Ymodel(j,:) = Model(mpar,pH);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
