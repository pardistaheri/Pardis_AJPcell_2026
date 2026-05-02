function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux_data = Data(:,1)/10^3; %Flux(mM/min)
A = Data(:,2)/10^3; %C_FUMm (M)
ATP = Data(:,3)/10^6; %C_ATPmm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux_data)
    Conc = [A(j), ATP(j)];
    Ymodel(j,:) = Flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux_data)./(Flux_data)).^2);
end
%Err = sum(((Cdata-Cmodel)./Cdata).^2);
