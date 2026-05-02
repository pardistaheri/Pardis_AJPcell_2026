clear all
close all
clc
%Driver
Data = load('Data.txt'); 
Flux_Data = Data(:,1); 
A_Data = Data(:,2)/10^3; %C_ACOAm (\muM)
B_Data = Data(:,3)/10^3; %C_OXAm (\muM)

mpar = load('mpar.txt');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
CIT_tot=[0:0.01:3];
for j=1:1:length(CIT_tot)
    A_tot = CIT_tot(j)/10^3;
    B_tot = 0;
    Conc=[A_tot, B_tot];
    flux1(:,j)=flux(mpar,Conc);
end

h1 = figure(1)
set(h1,'Position',[10 15 600 300])
subplot(1,2,1);
plot(A_Data(1:9,:),Flux_Data(1:9,:),'k*',CIT_tot,flux1(:,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('CIT (mM)')
ylabel('ACON Forward Flux (\muM/sec)')
text(0.2,0.155,'A','Fontsize',18)
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
ICIT_tot=[0:0.01:1.5];
for j=1:1:length(ICIT_tot)
    A_tot = 0;
    B_tot = ICIT_tot(j)/10^3;
    Conc=[A_tot, B_tot];
    flux2(:,j)=-flux(mpar,Conc);
end

subplot(1,2,2);
plot(B_Data(10:18,:),-Flux_Data(10:18,:),'k*',ICIT_tot,10^9*0.42*flux2(:,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ICIT (mM)')
ylabel('ACON reverse Flux (\muM/sec)')
text(0.1,0.3,'B','Fontsize',18)
box off