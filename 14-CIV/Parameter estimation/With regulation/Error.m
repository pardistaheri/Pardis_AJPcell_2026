function [Err] = Error(mpar)

% Load experimental data
Data = load('Data_C.txt');
Flux = Data(:,1)/10^3; %Flux(mmol/min)
fract = Data(:,2); 
B = Data(:,3)/10^6; %O2 (M)
dPsi = Data(:,4); %membrane potential (M)
pHi = Data(:,5); %inner-membrane pH (M)
dpH = Data(:,6); %difference of matrix and inner-membrane pH (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [fract(j), B(j)];
    Ymodel(j,:) = flux(mpar,dPsi(j),Conc,pHi(j),dpH(j));
end
%Ymodel=Rflux/60*1e6;
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end
