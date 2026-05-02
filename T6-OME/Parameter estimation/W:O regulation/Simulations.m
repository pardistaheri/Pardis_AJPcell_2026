%%%Load experimental data
Data = load('Data.txt'); 
Flux_Data = Data(:,1); %Flux(nM/min/ug)
A_data = Data(:,2); %C_AKGe (M)
B_data = Data(:,3); %C_MALm (M)
D_data = Data(:,5); %C_MALe (M)

%%%Load Parameters
mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
A_tot = [0:0.001:0.4]; 
for j=1:1:length(A_tot)
    A = A_tot(j)/10^3;
    B = 20/10^3;
    Conc=[A,B,0,0];
    flux1(j,:)=flux(mpar,Conc)*10^6;
end
h1 = figure(1);
set(h1,'Position',[10 10 900 300])
subplot(1,3,1)
plot(A_data(1:7,:),Flux_Data(1:7,:),'k*',A_tot,flux1,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('AKG_e (mM)')
ylabel('Influx of AKG(\mumol/mg/min)')
legend('MAL_m=20mM','Location','southeast')
%ylim([0 0.3])
text(0.01,6,'A','Fontsize',18)
box off
legend box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2 
B_tot = [0:0.1:4];
for i=1:1:length(B_tot)
        A = 0.5/10^3;
        B = B_tot(i)/10^3;
        Conc=[A, B, 0,0];
        flux2(i,:)=flux(mpar,Conc)*10^6;
end
subplot(1,3,2)
plot(B_data(8:13,:),Flux_Data(8:13,:),'k*',B_tot,1.3*flux2,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('MAL_m (mM)')
ylabel('Eflux of MAL(\mumol/mg/min)')
legend('AKG_e=0.5mM','Location','southeast')
text(0.1,6,'B','Fontsize',18)
box off
legend box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 3 
D_tot = [0:0.01:1.5];
for i=1:1:length(D_tot)
        A = 0;
        B = 0;
        C = 20/10^3;
        D = D_tot(i)/10^3;
        Conc=[A, B, C, D];
        flux3(i,:)=-flux(mpar,Conc)*10^6;
end
subplot(1,3,3)
plot(D_data(14:20,:),-Flux_Data(14:20,:),'k*',D_tot,flux3,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('MAL_e (mM)')
ylabel('Influx of MAL(\mumol/mg/min)')
legend('AKG_m=20mM','Location','southeast')
text(0.05,5,'C','Fontsize',18)
box off
legend box off
