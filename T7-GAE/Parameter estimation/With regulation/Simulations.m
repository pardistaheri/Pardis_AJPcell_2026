%%%Load experimental data
Data = load('Data.txt'); 
Flux_Data = Data(:,1); %Flux(nM/min/ug)
A_data = Data(:,2); %C_ASPm (M)

%%%Load Parameters
mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
A_tot = [0:0.01:10]; %ASPm
Ca_tot = [10 36 68 91 111];
for i=1:1:length(Ca_tot)
    for j=1:1:length(A_tot)
    A = A_tot(j)/10^3;
    B = 10/10^3;
    Ca = Ca_tot(i)/10^6;
   
    Conc=[A, B, Ca];
    flux1(i,j)=1.1*flux(mpar,Conc)*10^6;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(A_data(1:7,:),Flux_Data(1:7,:),'k*',A_data(8:16,:),Flux_Data(8:16,:),'ko',A_data(17:23,:),Flux_Data(17:23,:),'k<',A_data(24:29,:),Flux_Data(24:29,:),'kd',A_data(30:36,:),Flux_Data(30:36,:),'kp',A_tot,flux1(1,:),'k',A_tot,flux1(2,:),'k',A_tot,1.2*flux1(3,:),'k',A_tot,flux1(4,:),'k',A_tot,0.75*flux1(5,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ASP_m (mM)')
%ylabel('Eflux of ASP(nmol/mg)')
legend('Ca=0.01mM','Ca=0.04mM','Ca=0.07mM','Ca=0.09mM','Ca=0.1mM','location','northwest')
text(0.3,0.25,'C','Fontsize',18)
box off
legend box off
