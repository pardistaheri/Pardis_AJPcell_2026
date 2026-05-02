%error function for fmincon
function Err = Error(mpar)
%load data
Data = load('Data.txt');
Flux = Data(44:71,1)/10^6; %Flux(mmol/min)
A = Data(44:71,2)/10^3; %C_ICITm (M)
B = Data(44:71,3)/10^3; %C_NADm (M)
D = Data(44:71,4)/10^6; %C_NADHm (M) 
ADP = Data(44:71,5)/10^3; %C_ADPm (M)
ATP = Data(44:71,5)/10^3; %C_ADPm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j), D(j), ADP(j), ATP(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end