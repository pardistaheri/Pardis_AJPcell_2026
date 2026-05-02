
%%% load parameter
mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
Data = load('Data_A.txt'); 
Flux_Data = Data(:,1); %Flux(uM/min/ug)
B = Data(:,3); %C_O2 (M)

O2_tot=[0:0.01:16];
for i=1:1:length(O2_tot)
    B_tot = O2_tot(i)/10^6;
    Conc=[0.1, B_tot];
    dPsi =180;
    pHi = 7.35;
    dpH = 0.25;
    flux1(i,:)=flux(mpar,dPsi,Conc,pHi,dpH)*10^3;
end
h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1);
plot(B(1:28,:),Flux_Data(1:28,:),'k*',O2_tot,flux1,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('O_2 (\muM)')
ylabel('CIV Forward Flux (\mumol/min)')
text(0.5,0.8,'A','Fontsize',18)
ylim([0 0.8])
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2 
Data = load('Data_B.txt'); 
Flux_Data = Data(:,1); %Flux(1/s)
A = Data(:,2); %C_UQ (M)

dPsi_tot = [126 160 178];
dpH_tot = [0 0.25 0.4];
fract_tot=[0:0.001:0.8];

for i=1:1:length(dPsi_tot)
    for j=1:1:length(fract_tot)
        fract = fract_tot(j);
        B_tot = 200/10^6;
        dPsi = dPsi_tot(i);
        dpH = dpH_tot(i);
        pHi = 7.35;
        Conc=[fract, B_tot];
        flux2(i,j)=flux(mpar,dPsi,Conc,pHi,dpH)*10^3;
    end
end
subplot(1,2,2);
plot(A(1:6,:),Flux_Data(1:6,:),'k*',A(7:16,:),Flux_Data(7:16,:),'ko',A(17:30,:),Flux_Data(17:30,:),'k<',fract_tot(1:170),flux2(1,1:170),'k',fract_tot(1:550),flux2(2,1:550)/1.8,'k',fract_tot,flux2(3,:)/3,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('Fractional Cytochrome reduced')
ylabel('CIV Forward Flux (1/sec)')
legend('deltaPsi=126mV, deltapH=0','deltaPsi=160mV, deltapH=0.25','deltaPsi=178mV, deltapH=0.4','Location','northwest')
ylim([0 60]);
text(0.02,60,'B','Fontsize',18)
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


