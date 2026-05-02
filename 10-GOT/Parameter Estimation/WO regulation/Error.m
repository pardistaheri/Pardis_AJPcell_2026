%error function for fmincon
function Err = Error(mpar)
%load data
Data = load('Data.txt');
Flux = Data(:,1); %Flux(mmol/min)
A = Data(:,2)/10^3; %C_ASPm (M)
B = Data(:,3)/10^3; %C_AKGm (M)
C = Data(:,4)/10^3; %C_GLUm (M) 
D = Data(:,5)/10^3; %C_OXAm (M)

%calculate the model flux based on free concentration
for j=1:length(Flux)
    Conc = [A(j), B(j), C(j), D(j)];
    Ymodel(j,:) = flux(mpar,Conc);
end
Err = sum(((Ymodel-Flux)./(Flux)).^2);
end