function [Err] = Error_R(mpar)

% Load experimental data
Data = load('Data_R.txt');
Flux = Data(:,1)/10^9; %Flux(mM/min)
A = 0; %C_Sucm (M)
B = 0; %C_FADmm (M)
C = 10/10^3; %C_Fumm (M)
D = Data(:,3)/10^6; %C_FADH2mm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A, B, C, D(j)];
    Ymodel(j,:) = flux_R(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
%Err = sum(((Cdata-Cmodel)./Cdata).^2);
