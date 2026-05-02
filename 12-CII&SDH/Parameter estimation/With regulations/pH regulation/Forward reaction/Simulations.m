% 1. CI biochemical reaction model fit
% the flux for CITS reaction in this model was devolped in Audi & Dash lab
% previously(Xiao et al 2018)
% Model fitted to data from Bazil et al. 2010
% The units are as follows Concentration M, Flux mmol/nin, pH= 8.1, Volume =1 ml, Mass in Microg,5 parameters

clear all
close all
clc
%Driver
mpar = load('mpar.txt');
%figure 1 pH regulation
Data = load('Data_f.txt'); 
Flux_Data = Data(:,1); %Flux(mM/min/ug)
pH_data = Data(:,2); %pHm (M)
pH=[5:0.01:9];
for i=1:1:length(pH)
    pH_tot = pH(i);
    Conc=[pH_tot];
    flux1(i,:)=flux(mpar,Conc);
end
save('v_f.txt','flux1','-ascii');
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(pH,flux1,'k',pH_data(1:7,:),Flux_Data(1:7,:),'kd','LineWidth',2.5,'MarkerSize',10.0);
set(gcf,'color','w')
set(gca,'FontSize',12)
xlabel('pH')
ylabel('Vmaxf_P_r_i_m_e/ Vmaxf')
text((0.6)*0.05,(0.08)*0.9,'B','Fontsize',14)
box off

