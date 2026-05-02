function [Err] = Error(mpar)

% Load experimental data
Data = load('Data.txt');
Flux = Data(:,1); %Flux(1/s)
A = Data(:,2)/10^6; %ADP_e
B = 5/10^3; %ATP_m
C = 5/10^3; %ADP_m
D = Data(:,3)/10^6; %ATP_e
dPsi = Data(:,4); %difference of matrix and inner-membrane pH (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j),B,C,D(j),dPsi(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
