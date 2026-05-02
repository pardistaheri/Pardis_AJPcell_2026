%error function for fmincon
function Err = Error(mpar)
%load data
Data = load('Data.txt');
Flux = Data(1:48,1)/1000; %Flux(mM/min)
A = Data(1:48,2)/1000; %C_PYRm (M)
B = Data(1:48,3)/1000; %C_COAm (M)
C = Data(1:48,4)/1000; %C_NADm (M) 

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j), C(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end