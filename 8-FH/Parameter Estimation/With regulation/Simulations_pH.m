%Simulations
clear all
close all
clc

mpar = load('mpar_pH.txt');

Data = load('Data_pH.txt'); 
Flux_Data = Data(:,1); %Flux(arbitrary unit)
pH_data = Data(:,2); %pHm (M)

pH_tot=[4.5:0.01:8.5];
for i=1:1:length(pH_tot)
    pH = pH_tot(i);
    flux1(i,:)=3.2*Flux_pH(mpar,pH);
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(pH_tot,flux1,'k',pH_data(:,:),Flux_Data(:,:),'kd','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
ylabel('FH Forward Flux (normalized)')
box off

