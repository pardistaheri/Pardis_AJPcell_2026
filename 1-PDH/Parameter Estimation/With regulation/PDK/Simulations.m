close all
clear all
clc

Data = load('Data.txt');
Flux_data = Data(:,1);
A_data  = Data(:,5); %C_PYRm (mM)
ATP_data = Data(:,6); %C_ATPm (mM)
ADP_data = Data(:,7); %C_ADPm (mM)
pH_data = Data(:,2); %pH

mpar = load('mpar.txt');

pH_m=[5.0:0.01:8.5];
for i=1:1:length(pH_m)
    %substrate concentrations
    A = 2/1000; %C_PYRm (M)
    ATP = 0.1/1000;
    ADP = 0;
    pH = pH_m(i);                 % pH = 7.5 for Tsai et al. 1973
    Conc = [A, ATP, ADP,pH];
    Pflux(i,:) = Flux(mpar,Conc); 
end

h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1);
plot(pH_m,Pflux,'k',pH_data(1:9,:),Flux_data(1:9,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
ylabel('PDK activity (\DeltaA_3_4_0/min)')
text(5.25,0.35,'G','Fontsize',18)
box off



ATP_m=[0.02:0.01:0.2];
for i=1:1:length(ATP_m)
    %substrate concentrations
    A = 2/1000; %C_PYRm (M)
    ATP = ATP_m(i)/1000;
    ADP = 0;
    pH=7.5;
    Conc = [A, ATP, ADP,pH];
    Aflux1(i,:) = Flux(mpar,Conc);
end

for i=1:1:length(ATP_m)
    %substrate concentrations
    A_tot = 2/1000; %C_PYRm (M)
    ATP = ATP_m(i)/1000;
    ADP = 0.5/1000;
    pH=7.5;
    Conc = [A, ATP, ADP,pH];
    Aflux2(i,:) = Flux(mpar,Conc); 
end

for i=1:1:length(ATP_m)
    %substrate concentrations
    A_tot = 0.5/1000; %C_PYRm (M)
    ATP = ATP_m(i)/1000;
    ADP = 0;
    pH=7.5;
    Conc = [A, ATP, ADP,pH];
    Aflux3(i,:) = Flux(mpar,Conc); 
end
subplot(1,2,2);
plot(ATP_data(10:14,:),Flux_data(10:14,:),'k*',ATP_data(15:19,:),Flux_data(15:19,:),'ko',ATP_data(20:24,:),Flux_data(20:24,:),'k<',ATP_m,1.6*Aflux1,'k',ATP_m,1.05*Aflux2,'k',ATP_m,0.75*Aflux3,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ATP (mM)')
% ylabel('PDK activity (\DeltaA_3_4_0/min)')
legend('ATP','ATP+0.5mM ADP','ATP+0.5mM PYR','Location','southeast');
text(0.015,0.3,'H','Fontsize',18)
box off
legend box off

