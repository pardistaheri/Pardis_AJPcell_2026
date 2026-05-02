%error function for fmincon
function Err = Error(mpar)
%load data
Data = load('Data.txt');
Flux = Data(:,1)/10^3; %Flux(mmol/min)
A = Data(:,2)/10^3; %C_AKGm (M)
B = Data(:,3)/10^3; %C_COAm (M)
C = Data(:,4)/10^3; %C_NADm (M) 
D = Data(:,5)/10^6; %C_SCOAm (M)
E = Data(:,6)/10^3; %C_NADHm (M) 
%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j), C(j), D(j), E(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end