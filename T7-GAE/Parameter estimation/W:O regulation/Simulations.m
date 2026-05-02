%%%Load experimental data
Data = load('Data.txt'); 
Flux_Data = Data(:,1); %Flux(nM/min/ug)
A_data = Data(:,2); %C_ASPm (M)

%%%Load Parameters
mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
A_tot = [0:0.01:15]; 
B_tot = [1.5 2.5 5 10 15];
for i=1:1:length(B_tot)
    for j=1:1:length(A_tot)
    A = A_tot(j)/10^3;
    B = B_tot(i)/10^3;
   
    Conc=[A, B,0];
    flux1(i,j)=flux(mpar,Conc)*10^6;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1)
plot(A_data(1:11,:),Flux_Data(1:11,:),'k*',A_data(12:22,:),Flux_Data(12:22,:),'ko',A_data(23:33,:),Flux_Data(23:33,:),'k<',A_data(34:43,:),Flux_Data(34:43,:),'kd',A_data(44:53,:),Flux_Data(44:53,:),'kp',...
    A_tot,flux1(1,:),'k',A_tot,flux1(2,:),'k',A_tot,1.05*flux1(3,:),'k',A_tot,flux1(4,:),'k',A_tot,flux1(5,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ASP_m (mM)')
ylabel('Eflux of ASP(nmol/mg)')
legend('GLU_e=1.5mM','GLU_e=2.5mM','GLU_e=5mM','GLU_e=10mM','GLU_e=15mM','Location','northwest')
ylim([0 0.3])
text(0.5,0.3,'A','Fontsize',18)
box off
legend box off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2 
Data = load('Data_A.txt'); 
Flux_Data = Data(:,1); %Flux(nM/min/ug)
dPsi_data = Data(:,4); %dPsi (M)

dPsi_tot = [10:0.1:200]; %Cae
for i=1:1:length(dPsi_tot)
        A = 10/10^3;
        B = 10/10^3;
        dPsi = dPsi_tot(i);
        Conc=[A, B, dPsi];
        flux2(i,:)=flux(mpar,Conc)*10^6;
end
subplot(1,2,2)
plot(dPsi_tot,flux2/20,'k',dPsi_data(1:13,:),Flux_Data(1:13,:)/60,'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('dPsi')
%ylabel('Eflux of ASP(nmol/mg)')
ylim([0 0.4])
text(5,0.4,'B','Fontsize',18)
box off
