%error function for fmincon
function Err = Error(mpar)
%load data
Data = load('Data.txt');
Flux = Data(:,1); %Flux(mM/min)
A = Data(:,2)/10^6; %C_ACOAm (M)
B = Data(:,3)/10^6; %C_OXAm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end