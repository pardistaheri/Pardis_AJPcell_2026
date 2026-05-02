clear all
close all
clc
%Driver
%load data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(U/ml)
A_data = Data(:,2); %C_GDPm (M)
B_data = Data(:,3); %C_ATPm (M)

mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
GDP_tot=[0:0.001:1.5];
    for j=1:1:length(GDP_tot)
    A_tot = GDP_tot(j)/10^3; %C_GTP (M)
    B_tot = 2/10^3; %C_ADP (M)
%     C_tot = ; %C_GDP (M)
%     D_tot = ; %C_ATP (M)
    Conc=[A_tot, B_tot, 0,0];
    flux1(:,j)=flux(mpar,Conc);
    end
    
h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1)
plot(A_data(1:12,:),Flux_Data(1:12,:),'k*',GDP_tot,0.9*flux1,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('GDP (mM)')
ylabel('NDK Reverse Flux (U/ml)')
text(0.08,0.01,'A','Fontsize',18)
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
ATP_tot=[0:0.001:4.5];
    for j=1:1:length(ATP_tot)
    A_tot = 0.7/10^3 ; %C_GTP (M)
    B_tot = ATP_tot(j)/10^3; %C_ADP (M)
%     C_tot = ; %C_GDP (M)
%     D_tot = ; %C_ATP (M)
    Conc=[A_tot, B_tot, 0,0];
    flux2(:,j)=flux(mpar,Conc);
    end
subplot(1,2,2)
plot(B_data(13:22,:),Flux_Data(13:22,:),'k*',ATP_tot,2.1*flux2,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
text(0.2,0.025,'B','Fontsize',18)
xlabel('ATP (mM)')
box off
