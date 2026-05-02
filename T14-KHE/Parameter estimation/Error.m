function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(1:15,1)/10^6; %Flux(mmol/min)
A = 0; %Km
B = 10.^-Data(1:15,2); %He
C = Data(1:15,3)/10^3; %Ke
D = Data(1:15,4)/10^9; %Hm

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A, B(j), C(j), D(j)];
    Ymodel(j,:) = Model(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
