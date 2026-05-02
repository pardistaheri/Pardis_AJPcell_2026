function [Err] = Error(mpar)

% Load experimental data
Data = load('Data_B.txt');
Flux = Data(:,1)/10^6; %Flux(mmol/min)
B = Data(:,2)/10^9;%He
C = Data(:,3)/10^3;%Nae
D = Data(:,4)/10^9;%Hm

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [B(j), C(j), D(j)];
    Ymodel(j,:) = Model(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
