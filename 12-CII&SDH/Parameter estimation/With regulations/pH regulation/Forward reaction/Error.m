function [Err] = Error(mpar)

% Load experimental data
Data = load('Data_f.txt');
Flux = Data(:,1); %Flux(mM/min)
pH = Data(:,2); %C_NADHm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    %Conc = [pH(j)];
    Ymodel(j,:) = flux(mpar,pH(j));
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
