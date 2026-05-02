function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(20:39,1)/10^3; %Flux(mM/min)
A = Data(20:39,2)/10^3; %C_Sucm (M)
B = Data(20:39,3)/10^3; %C_FADmm (M)
C = Data(20:39,4)/10^3; %C_Fumm (M)
D = Data(20:39,5)/10^3; %C_FADH2mm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j), C(j), D(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
%Err = sum(((Cdata-Cmodel)./Cdata).^2);
