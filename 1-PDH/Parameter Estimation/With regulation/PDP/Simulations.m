close all
clear all
clc

Data = load('Data.txt');
Flux_data = Data(:,1);
pH_data = Data(:,2);
Mg_data = Data(:,3);          
Ca_data = Data(:,4)/1000;         
mpar = load('mpar.txt');

pH_m=[5.0:0.01:8.5];
for i=1:1:length(pH_m)
    %%% Experimental conditions (total cations)
    pH_tot = pH_m(i);                 
    H_tot = 10.^-pH_tot;                
    Mg_tot = 10/1000;          
    Ca_tot = 0;          
    Na_tot = 0;                     
   
    Conc = [H_tot, Mg_tot, Ca_tot];
    Phflux(i,:) = Flux(Conc,mpar);
    Pflux(i,:)=Phflux(i,:);
end

subplot(1,3,1);
plot(pH_m,1.15*Pflux,'k',pH_data(1:9,:),Flux_data(1:9,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
ylabel('PDP activity (\DeltaA_3_4_0/min)')
text(5.25,0.5,'D','Fontsize',18)
box off

Mg_m=[0:0.1:10];
for i=1:1:length(Mg_m)
    % %% Experimental conditions (total cations)
    pH_tot = 7;                 
    H_tot = 10.^-pH_tot;                 
    Mg_tot = Mg_m(i)/1000;          
    Ca_tot = 0;          
    Na_tot = 0;                       
  
    Conc = [H_tot, Mg_tot, Ca_tot];
    Mflux(i,:) = Flux(Conc,mpar); 
end

h1 = figure(1)
set(h1,'Position',[10 10 900 300])
subplot(1,3,2);
plot(Mg_m,Mflux,'k',Mg_data(10:21,:),Flux_data(10:21,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
text(0.62,0.5,'E','Fontsize',18)
xlabel('Mg (mM)')
box off

Ca_m=[0:0.001:0.6];
for i=1:1:length(Ca_m)
    % %% Experimental conditions (total cations)
    pH_tot = 7;                 
    H_tot = 10.^-pH_tot;               
    Mg_tot = 0/1000;        
    Ca_tot = Ca_m(i)/1000;         
    Na_tot = 0;                     
    
    Conc = [H_tot, Mg_tot, Ca_tot];
    Cflux(i,:) = Flux(Conc,mpar);
end
subplot(1,3,3);
plot(Ca_m,Cflux,'k',Ca_data(22:27,:),Flux_data(22:27,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('Ca (mM)')
text(0.05,0.12,'F','Fontsize',18)
ylim([0 0.12])
box off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
