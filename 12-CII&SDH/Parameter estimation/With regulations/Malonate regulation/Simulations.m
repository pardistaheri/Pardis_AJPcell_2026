% Load experimental data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(mM/min)
A = Data(:,2); %C_Sucm (uM)

mpar = load('mpar.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
SUC_tot=[0:0.01:5]; %
MA=[0 8.9 17.9]*10^-6; %
for i=1:1:length(MA)
    for j=1:1:length(SUC_tot)
        A_tot = SUC_tot(j)*10^-3;
        MA_tot = MA(i);
    Conc=[A_tot, MA_tot];
    flux1(i,j)=flux(mpar,Conc);
   
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(A(1:4,:),Flux_Data(1:4,:),'k*',A(5:8,:),Flux_Data(5:8,:),'ko',A(9:11,:),Flux_Data(9:11,:),'kd',...
    SUC_tot,flux1(1,:),'k',SUC_tot,flux1(2,:),'k',SUC_tot,flux1(3,:),'k', 'LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('SUC (mM)')
ylabel('CII/SDH Forward Flux (mmol/min)')
legend('0 \muM Malonate', '8.9 \muM Malonate','17.9 \muM Malonate','Location','southeast')
text(0.4,0.06,'D','Fontsize',18)
legend boxoff
box off

