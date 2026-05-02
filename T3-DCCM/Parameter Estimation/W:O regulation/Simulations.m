%%%Load experimental data
Data = load('Data.txt'); 
Flux_Data = Data(:,1); %Flux(mM/min)
B_data = Data(:,3); %C_MALe (M)
C_data = Data(:,4); %C_Pie (M)

%%%Load Parameters
mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
C_tot = [0:0.01:1.05]; %Cae
for i=1:1:length(C_tot)
    A = 0;
    B = 0;
    C = C_tot(i)/10^3;
    D = 15/10^3;
   
    Conc=[A, B, C, D];
    flux1(:,i)=-flux(mpar,Conc);
end

% h1 = figure(1);
% set(h1,'Position',[10 10 900 300])
% subplot(1,3,1)

h1 = figure(1);
set(h1,'Position',[10 10 300 300])
plot(C_data(5:8,:),-Flux_Data(5:8,:),'k*',C_tot,2.25*flux1,'k','LineWidth',1.5,'MarkerSize',14.0);
set(gcf,'color','w')
%set(gca,'FontSize',14,'LineWidth',2.0)
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',2.0)
xlabel('Pi_e (mM)','FontSize',14)
ylabel('Eflux of Malate (mmol/min/mg)','FontSize',14)
legend('MAL_m = 15mM','Location','southeast')
text(0.05,2.5,'A','FontWeight','bold','Fontsize',18)
box off
legend box off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2 
B_tot = [0:0.001:1]; %Cae
for i=1:1:length(B_tot)
    A = 15/10^3;
    B = B_tot(i)/10^3;
    C = 0;
    D = 0;
   
    Conc=[A, B, C, D];
    flux2(:,i)=flux(mpar,Conc);
end

% subplot(1,3,2)
h2 = figure(2);
set(h2,'Position',[10 10 300 300])
plot(B_data(1:4,:),Flux_Data(1:4,:),'k*',B_tot,0.9*flux2,'k','LineWidth',1.5,'MarkerSize',14.0);
set(gcf,'color','w')
%set(gca,'FontSize',14,'LineWidth',2.0)
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',2.0)
xlabel('MAL_e (mM)','FontSize',14)
ylabel('Influx of Malate (mmol/min/mg)','FontSize',14)
legend('Pi_m = 15mM','Location','southeast')
text(0.05,4,'B','FontWeight','bold','Fontsize',18)
box off
legend box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 3 
Data = load('Data_A.txt'); 
Flux_Data = Data(:,1); %Flux(mM/min)
B_data = Data(:,3); %C_MALe (M)

C_tot = [4 3 2 1 0];
B_tot = [0:0.001:1.0]; %Cae
for i=1:1:length(C_tot)
    for j=1:1:length(B_tot)
        A = 15/10^3;
        B = B_tot(j)/10^3;
        C = C_tot(i)/10^3;
        D = 0;
   
    Conc=[A, B, C, D];
    flux3(i,j)=flux(mpar,Conc);
    end
end
% subplot(1,3,3)
h3 = figure(3);
set(h3,'Position',[10 10 300 300])
plot(B_data(1:4,:),Flux_Data(1:4,:),'k*',B_data(5:8,:),Flux_Data(5:8,:),'ko',B_data(9:12,:),Flux_Data(9:12,:),'k<',B_data(13:16,:),Flux_Data(13:16,:),'kd',B_data(17:20,:),Flux_Data(17:20,:),'kp',...
    B_tot,flux3(1,:),'k',B_tot,flux3(2,:),'k',B_tot,flux3(3,:),'k',B_tot,flux3(4,:),'k',B_tot,flux3(5,:),'k','LineWidth',1.5,'MarkerSize',14.0);
set(gcf,'color','w')
%set(gca,'FontSize',14,'LineWidth',2.0)
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',2.0)
xlabel('MAL_e (mM)','FontSize',14)
ylabel('Influx of Malate (mmol/min/mg)','FontSize',14)
legend('Pi_e = 4mM','Pi_e = 3mM','Pi_e = 2mM','Pi_e = 1mM','Pi_e = 0mM','Location','northwest')
ylim([0 6])
text(0.05,6,'C','FontWeight','bold','Fontsize',18)
box off
legend box off