%Simulation for reverse reaction
%Driver
clear all
close all 
clc

mpar = load('mpar.txt');
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
Data = load('Data_Rev.txt'); 
Flux_Data = Data(:,1); %Flux(muM/min)
C = Data(:,4); %C_OXAm (M)

OXA_tot=[0:0.001:0.055]; %
NADH_tot=[0.0067 0.01 0.015 0.02 0.025 0.05]; %

for i=1:1:length(NADH_tot)
    for j=1:1:length(OXA_tot)
        A_tot = 0;
        B_tot = 0;
        C_tot = OXA_tot(j)/10^3;
        D_tot = NADH_tot(i)/10^3;
    Conc=[A_tot, B_tot,C_tot, D_tot];
    flux1(i,j)=-0.6*Model(mpar,Conc)*10^3;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 400 400])
plot(OXA_tot,flux1(1,:),'g',OXA_tot,flux1(2,:),'m',OXA_tot,flux1(3,:),'c',OXA_tot,flux1(4,:),'k',OXA_tot,1.1*flux1(5,:),'r',OXA_tot,1.25*flux1(6,:),'b',...
    C(1:7,:),Flux_Data(1:7,:),'g*',C(8:14,:),Flux_Data(8:14,:),'m*',C(15:21,:),Flux_Data(15:21,:),'c*',C(22:28,:),Flux_Data(22:28,:),'k*',C(29:35,:),Flux_Data(29:35,:),'r*',C(36:42,:),Flux_Data(36:42,:),'b*','LineWidth',1,'MarkerSize',5);set(gcf,'color','w')
set(gca,'FontSize',10)
xlabel('OXA(mM)')
ylabel('Reverse Flux (\mumol/min/\mug)')
legend('0.0067 mM NADH','0.01 mM NADH','0.015 mM NADH','0.02 mM NADH','0.025 mM NADH','0.05 mM NADH','Location','northwest')
text('Fontsize',16)
title('MAL=0mM,NAD=0mM')
legend boxoff
box off