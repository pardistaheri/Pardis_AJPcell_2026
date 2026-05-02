%error function for fmincon
function Err = Error(mpar)
Data = load('Data.txt');
Flux_data = Data(:,1);
A  = Data(:,5)/1000; %C_PYRm (M)
ATP = Data(:,6)/1000;
ADP = Data(:,7)/1000;
pH = Data(:,2);

for j=1:1:length(Flux_data)
    Conc = [A(j), ATP(j), ADP(j),pH(j)];
    Ymodel(j,:) = Flux(mpar,Conc);  
end
%calculate the model flux based on free concentration
    Err = sum(((Ymodel-Flux_data)./(Flux_data)).^2);
end
