clear all
close all
clc
%Driver
mpar = load('mpar_pH.txt');
%figure 1 pH regulation
Data = load('Data_pH.txt'); 
Flux_Data = Data(:,1); %Flux(mM/min/ug)
pH_data = Data(:,2); %pHm (M)
pH=[5:0.01:9.5];
for i=1:1:length(pH)
    pH_tot = pH(i);
    Conc=[pH_tot];
    flux1(i,:)=flux_pH(mpar,Conc);
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(pH,flux1,'k',pH_data(1:18,:),Flux_Data(1:18,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
ylabel('Vmaxf (normalized)')
text(5.2,1,'E','Fontsize',18)
xlim([5 9.5])
box off

