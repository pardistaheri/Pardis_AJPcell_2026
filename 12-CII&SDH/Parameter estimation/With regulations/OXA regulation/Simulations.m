% Load experimental data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(mM/min)
A = Data(:,2); %C_Sucm (uM)

mpar = load('mpar.txt');
options = [];       % ODE options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
SUC_tot=[0:0.01:3.5]; %
OAA=[0 5.4 10.8]*10^-6; %
for i=1:1:length(OAA)
    for j=1:1:length(SUC_tot)
        A_tot = SUC_tot(j)*10^-3;
        OAA_tot = OAA(i);
    Conc=[A_tot, OAA_tot];
    flux1(i,j)=flux(mpar,Conc)*10^3;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(A(1:5,:),Flux_Data(1:5,:),'k*',A(6:10,:),Flux_Data(6:10,:),'ko',A(11:15,:),Flux_Data(11:15,:),'kd',...
    SUC_tot,flux1(1,:),'k',SUC_tot,flux1(2,:),'k',SUC_tot,flux1(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('SUC (mM)')
ylabel('CII/SDH Forward Flux (\mumol/min)')
legend('0 \muM OAA', '5.4 \muM OAA','10.8 \muM OAA','Location','northwest')
text(0.2,2,'C','Fontsize',18)
%ylim([0 2.5])
legend boxoff
box off
