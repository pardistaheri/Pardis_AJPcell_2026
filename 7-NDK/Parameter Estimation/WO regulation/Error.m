%error function for fmincon
function Err = Error(mpar)
%load data
Data = load('Data.txt');
Flux = Data(:,1); %Flux(mmol/min)
A = Data(:,2)/10^3; %C_GTPm (M)
B = Data(:,3)/10^3; %C_ADPm (M)
C = 0; %C_GDPm (M) 
D = 0; %C_ATPm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j), C, D];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end