function [Err] = Error_pH(mpar)

% Load experimental data
Data = load('Data_pH.txt');
Flux_data = Data(:,1); %Flux(arbitrary unit)
pH_data = Data(:,2); 

%calculate the model flux based on free concentration
for j=1:length(Flux_data)
    Ymodel(j,:) = flux_pH(mpar,pH_data(j));
end
Err = sum(((Ymodel-Flux_data)./(Flux_data)).^2);
end
