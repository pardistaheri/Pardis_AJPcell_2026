clear all
close all
clc
%Driver
%load data
Data = load('Data.txt');
Flux_Data = Data(:,1)*10^3; %Flux(nmol/min)
A_data = Data(:,2); %C_ASPm (M)
B_data = Data(:,3); %C_AKGm (M)

mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
ASP_tot=[0:0.001:5];
AKG_tot=[0.0833 0.111 0.167 0.333];
for i=1:1:length(AKG_tot)
    for j=1:1:length(ASP_tot)
    A_tot = ASP_tot(j)/10^3; %C_ASP (M)
    B_tot = AKG_tot(i)/10^3; %C_AKG (M)
    C_tot = 0; %C_GLU (M)
    D_tot = 0; %C_OXA (M)
    Conc=[A_tot, B_tot, C_tot,D_tot];
    flux1(i,j)=flux(mpar,Conc)*10^3;
    end
end

h1 = figure(1)
set(h1,'Position',[10 10 600 600])
subplot(2,2,1)
plot(A_data(1:4,:),Flux_Data(1:4,:),'k*',A_data(5:8,:),Flux_Data(5:8,:),'k<',A_data(9:12,:),Flux_Data(9:12,:),'ko',A_data(13:16,:),Flux_Data(13:16,:),'kd',ASP_tot,1.8*flux1(1,:),'k',ASP_tot,1.8*flux1(2,:),'k',ASP_tot,1.8*flux1(3,:),'k',ASP_tot,1.8*flux1(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ASP (mM)')
ylabel('GOT Forward Flux (1/\mug/ml/min)')
legend('0.0833 mM AKG','0.111 mM AKG','0.167 mM AKG','0.333 mM AKG','Location','southeast')
text(0.2,6,'A','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
AKG_tot=[0:0.001:3];
GLU_tot=[0 8 16 24];
for i=1:1:length(GLU_tot)
    for j=1:1:length(AKG_tot)
    A_tot = 10/10^3; %C_ASP (M)
    B_tot = AKG_tot(j)/10^3; %C_AKG (M)
    C_tot = GLU_tot(i)/10^3; %C_GLU (M)
    D_tot = 0; %C_OXA (M)
    Conc=[A_tot, B_tot, C_tot,D_tot];
    flux3(i,j)=flux(mpar,Conc)*10^3;
    end
end
subplot(2,2,2)
plot(B_data(17:20,:),Flux_Data(17:20,:),'k*',B_data(21:24,:),Flux_Data(21:24,:),'k<',B_data(25:28,:),Flux_Data(25:28,:),'ko',B_data(29:32,:),Flux_Data(29:32,:),'kd',AKG_tot,0.8*flux3(1,:),'k',AKG_tot,0.95*flux3(2,:),'k',AKG_tot,flux3(3,:),'k',AKG_tot,flux3(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('AKG (mM)')
legend('0 mM GLU','8 mM GLU','16 mM GLU','24 mM GLU','Location','southeast')
text(0.1,5,'B','Fontsize',18)
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 3 
ASP_tot=[0:0.01:20];
OXA_tot=[0 0.08 0.16 0.24];
for i=1:1:length(OXA_tot)
    for j=1:1:length(ASP_tot)
    A_tot = ASP_tot(j)/10^3; %C_ASP (M)
    B_tot = 1.33/10^3; %C_AKG (M)
    C_tot = 0; %C_GLU (M)
    D_tot = OXA_tot(i)/10^3; %C_OXA (M)
    Conc=[A_tot, B_tot, C_tot,D_tot];
    flux4(i,j)=flux(mpar,Conc)*10^3;
    end
end
subplot(2,2,3)
plot(A_data(33:36,:),Flux_Data(33:36,:),'k*',A_data(37:40,:),Flux_Data(37:40,:),'k<',A_data(41:44,:),Flux_Data(41:44,:),'ko',A_data(45:48,:),Flux_Data(45:48,:),'kd',ASP_tot,0.65*flux4(1,:),'k',ASP_tot,0.8*flux4(2,:),'k',ASP_tot,flux4(3,:),'k',ASP_tot,flux4(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ASP (mM)')
ylabel('GOT Forward Flux (1/\mug/ml/min)')
legend('0 mM OXA','0.08 mM OXA','0.16 mM OXA','0.24 mM OXA','Location','southeast')
text(0.5,4,'C','Fontsize',18)
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 4 
AKG_tot=[0:0.001:3];
OXA_tot=[0 0.04 0.08 0.12];
for i=1:1:length(OXA_tot)
    for j=1:1:length(AKG_tot)
    A_tot = 40/10^3; %C_ASP (M)
    B_tot = AKG_tot(j)/10^3; %C_AKG (M)
    C_tot = 0; %C_GLU (M)
    D_tot = OXA_tot(i)/10^3; %C_OXA (M)
    Conc=[A_tot, B_tot, C_tot,D_tot];
    flux6(i,j)=1.2*flux(mpar,Conc)*10^3;
    end
end
subplot(2,2,4)
plot(B_data(49:52,:),Flux_Data(49:52,:),'k*',B_data(53:56,:),Flux_Data(53:56,:),'k<',B_data(57:60,:),Flux_Data(57:60,:),'ko',B_data(61:64,:),Flux_Data(61:64,:),'kd',AKG_tot,flux6(1,:),'k',AKG_tot,1.1*flux6(2,:),'k',AKG_tot,1.2*flux6(3,:),'k',AKG_tot,1.2*flux6(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('AKG (mM)')
legend('0 mM OXA','0.04 mM OXA','0.08 mM OXA','0.12 mM OXA','Location','southeast')
text(0.1,8,'D','Fontsize',18)
legend boxoff
box off
