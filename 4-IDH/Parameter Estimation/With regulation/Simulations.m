clear all
close all
clc
%Driver
%load data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(nmol/min)
A_data = Data(:,2); %C_ICITm (M)
B_data = Data(:,3); %C_NADm (M)
D_data = Data(:,4); %C_NADHm (M) 

mpar = load('mpar.txt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
ADP_tot = [0 1.5];
ATP_tot = [1.5 0];
ICIT_tot=[0:0.01:2.5];
for i=1:1:length(ADP_tot)
    for j=1:1:length (ICIT_tot)
    A_tot = ICIT_tot(j)/1000; %C_ICITm (M)
    B_tot = 2/10^3; %C_NADm (M)
    D_tot = 0; %C_NADHm (M)
    ADP = ADP_tot(i)/10^3;
    ATP = ATP_tot(i)/10^3;
    Conc=[A_tot, B_tot,D_tot,ADP,ATP];
    flux1(i,j)=10^6*flux(mpar,Conc);
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1);
plot(A_data(44:55,:),Flux_Data(44:55,:),'k*',A_data(56:71,:),Flux_Data(56:71,:),'ko',ICIT_tot,flux1(1,:),'k',ICIT_tot,flux1(2,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ICIT (mM)')
ylabel('IDH Forward Flux (munit/mg protein)')
legend('ATP=1.5 mM','ADP= 1.5 mM','Location','southeast')
text(0.2,60,'D','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 6 
mpar_pH = load('mpar_pH.txt');
%figure 1 pH regulation
Data_pH = load('Data_pH.txt'); 
Flux_Data = Data_pH(:,1); %Flux(arbitrary unit)
pH_data = Data_pH(:,2); %pHm (M)

pH_tot=[5:0.01:9];
for i=1:1:length(pH_tot)
    pH = pH_tot(i);
    flux2(i,:)=flux_pH(mpar_pH,pH);
end
subplot(1,2,2);
plot(pH_tot,flux2,'k',pH_data(:,:),Flux_Data(:,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
ylabel('IDH Forward Flux (arbitrary units)')
text(5.2,0.4,'E','Fontsize',18)
box off
