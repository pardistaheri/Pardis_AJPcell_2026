clear all
close all
clc
%Driver
%load data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(nmol/min)
A_data = Data(:,2); %C_ICITm (M)
Ca_data = Data(:,5); %C_NADm (M)

mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
Ca_tot = [0 91.21];
AKG_tot=[0:0.01:10];
for i=1:1:length(Ca_tot)
    for j=1:1:length (AKG_tot)
    A_tot = AKG_tot(j)/10^3; %C_AKGDHm (M)
    B_tot = 0.25/10^3; %C_COAm (M)
    C_tot = 0.1/10^3; %C_NADm (M)
    Ca = Ca_tot(i)/10^6;
    Conc=[A_tot, B_tot,C_tot,Ca];
    flux1(i,j)=flux(mpar,Conc);
    end
end
   
h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1);
plot(A_data(1:11,:),Flux_Data(1:11,:),'ko',A_data(12:22,:),Flux_Data(12:22,:),'k*',AKG_tot,2.8*flux1(1,:),'k',AKG_tot,2.8*flux1(2,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('AKG (mM)')
ylabel('AKGDH Forward Flux (unit/mg protein)')
legend('Ca^2^+ = 0mM','Ca^2^+ =90\muM','Location','northeast')
ylim([0 20])
text(0.5,20,'I','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
Ca_tot3 = [0.001:0.001:100];

for i=1:1:length(Ca_tot3)
    A_tot = 0.1/10^3; %C_AKGDHm (M)
    B_tot = 0.25/10^3; %C_COAm (M)
    C_tot = 1/10^3; %C_NADm (M)
    Ca = Ca_tot3(i)/10^6;
    Conc=[A_tot, B_tot,C_tot,Ca];
    flux3(i,:)=flux(mpar,Conc);
end
    
subplot(1,2,2);
plot(log(Ca_data(23:34,:)),Flux_Data(23:34,:),'k*',log(Ca_tot3),flux3(:,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('Ca^2^+ (\muM)')
ylabel('AKGDH Forward Flux (nmol/mg protein)')
legend('AKG = 0.1mM','Location','southeast')
text(-6.5,5,'J','Fontsize',18)
legend boxoff
box off
