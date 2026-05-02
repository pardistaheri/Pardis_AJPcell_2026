% 1. CI biochemical reaction model fit
% the flux for CITS reaction in this model was devolped in Audi & Dash lab
% previously(Xiao et al 2018)
% Model fitted to data from Bazil et al. 2010
% The units are as follows Concentration M, Flux mmol/nin, pH= 8.1, Volume =1 ml, Mass in Microg,5 parameters
close all
clear all
clc

%Driver
mpar = load('mpar.txt');
%figure 1 pH regulation
Data = load('Data_pH.txt'); 
Flux_Data = Data(:,1); %Flux%
pH_data = Data(:,2); %pHm 
pH=[6.9:0.01:7.5];
for i=1:1:length(pH)
    pH_tot = pH(i);
    flux1(i,:)=Model(mpar,pH_tot)*10^6;
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(pH,flux1,'k',pH_data(1:7,:),Flux_Data(1:7,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH_m')
ylabel('Eflux(nmol/s/mg)')
text(6.82,6,'C','Fontsize',16)
box off

