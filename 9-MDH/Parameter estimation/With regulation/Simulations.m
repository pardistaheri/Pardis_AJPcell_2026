% Load experimental data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(mM/min)
D = Data(:,2); %C_NADHm (mM)

mpar = load('mpar.txt');
options = [];       % ODE options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
NADH_tot=[0:0.0001:0.025]; 
for i=1:1:length(NADH_tot)
    D_tot = NADH_tot(i)*10^-3;
    ATP_tot = 2/10^3;
    ADP_tot = 0;
    AMP_tot = 0;
    Conc=[D_tot,ATP_tot,ADP_tot,AMP_tot];
    flux1(i,:)=flux(mpar,Conc)*10^3;
end
for i=1:1:length(NADH_tot)
    D_tot = NADH_tot(i)*10^-3;
    ATP_tot = 0;
    ADP_tot = 2/10^3;
    AMP_tot = 0;
    Conc=[D_tot,ATP_tot,ADP_tot,AMP_tot];
    flux2(i,:)=flux(mpar,Conc)*10^3;
end
for i=1:1:length(NADH_tot)
    D_tot = NADH_tot(i)*10^-3;
    ATP_tot = 0;
    ADP_tot = 0;
    AMP_tot = 2/10^3;
    Conc=[D_tot,ATP_tot,ADP_tot,AMP_tot];
    flux3(i,:)=flux(mpar,Conc)*10^3;
end
for i=1:1:length(NADH_tot)
    D_tot = NADH_tot(i)*10^-3;
    ATP_tot = 0;
    ADP_tot = 0;
    AMP_tot = 0;
    Conc=[D_tot,ATP_tot,ADP_tot,AMP_tot];
    flux4(i,:)=flux(mpar,Conc)*10^3;
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(D(1:4,:),Flux_Data(1:4,:),'k*',D(5:8,:),Flux_Data(5:8,:),'ko',D(9:13,:),Flux_Data(9:13,:),'kd',D(14:28,:),Flux_Data(14:28,:),'k<',...
    NADH_tot,flux1,'k',NADH_tot,flux2,'k',NADH_tot,flux3,'k',NADH_tot,flux4,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
legend('2 mM ATP','2 mM ADP','2 mM AMP','W/O AXP')
xlabel('NADH (mM)')
ylabel('MDH Reverse Flux (\mumol/min)')
text(0.002,0.2,'G','Fontsize',18)
ylim([0 0.2])
box off
legend box off

