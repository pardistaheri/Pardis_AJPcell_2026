% ComplexII/Succinate dehydregenase parameter estimation forward and
% reverse reaction with pH regulation
% SUC + UQ = FUM + UQH2
% A - SUC; B - UQ; C - FUM; D -UQH2 ;

% NOTE: The following code was adopted from the lung model developed
% by Zhang et al 2018 published in Public Library of Science 
% Data used for parameter estimations were obtained from Maklashina et al
% 1999
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024

clear all
close all
clc

Data_r = load('Data_r.txt'); 
Flux_r = Data_r(:,1); %Flux(mM/min/ug)
pH_r = Data_r(:,2); %pHm (M)

Data_f = load('Data_f.txt'); 
Flux_f = Data_f(:,1); %Flux(mM/min/ug)
pH_f = Data_f(:,2); %pHm (M)

v_f = load('v_r.txt'); 
V_r = v_f(:,1); %Flux(mM/min/ug)

v_f = load('v_f.txt'); 
V_f = v_f(:,1); %Flux(mM/min/ug)

pH=[5:0.01:9];
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(pH_f(:,:),Flux_f(:,:),'k*',pH_r(:,:),Flux_r(:,:),'ko',pH,V_f,'k',pH,V_r,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
ylabel('Vmax (normalized)')
ylim([0 1.2])
legend('V_f_o_r_w_a_r_d', 'V_r_e_v_e_r_s_e','Location','southeast')
text(5.2,1.2,'E','Fontsize',18)
legend box off
box off
