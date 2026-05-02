% 1. CI biochemical reaction model fit
% the flux for CITS reaction in this model was devolped in Audi & Dash lab
% previously(Xiao et al 2018)
% Model fitted to data from Bazil et al. 2010
% The units are as follows Concentration M, Flux mmol/nin, pH= 8.1, Volume =1 ml, Mass in Microg,5 parameters

clear all
close all
clc
%Driver
mpar = load('mpar_pH.txt');
%figure 1 pH regulation
Data = load('Data_pH.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux(mM/min/ug)
pH_data = Data(:,2); %pHm (M)
pH=[6:0.01:9];
for i=1:1:length(pH)
    pH_tot = pH(i);
    Conc=[pH_tot, 0];
    flux1(i,:)=flux_pH(mpar,Conc)*10^3;
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(pH,flux1,'k',pH_data(1:5,:),Flux_Data(1:5,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
ylabel('Flux (\mumol/min)')
ylim([0.2 1.2])
box off

