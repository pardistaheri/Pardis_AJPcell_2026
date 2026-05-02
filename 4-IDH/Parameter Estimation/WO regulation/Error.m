%error function for fmincon
function Err = Error(mpar)
%load data
Data = load('Data.txt');
Flux = Data(:,1)/10^6; %Flux(mmol/min)
A = Data(:,2)/10^3; %C_ICITm (M)
B = Data(:,3)/10^3; %C_NADm (M)
D = Data(:,4)/10^6; %C_NADHm (M) 

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j), D(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end