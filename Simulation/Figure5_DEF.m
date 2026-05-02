%% clear workspace 
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read param
run Read_Condt
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax
p.Ve=1.33e-3;

pyr_data = load('PYR.txt'); 
mal_data = load('MAL.txt'); 
cit_data = load('CIT.txt'); 
Time_Data = pyr_data(:,1); %Flux(mM/min/ug)
PYR_Data = pyr_data(:,2); %C_PYRm (M)
MAL_Data = mal_data(:,2); %C_PYRm (M)
CIT_Data = cit_data(:,2); %C_PYRm (M)
%% substrates and ADP additions 
PYR_index=[0 0.5 1 5 10]/10^3;
MAL_index=5/10^3;
p.NSub=length(PYR_index);

ADP_add=[0 0 0]*1e-6; p.ADPL=length(ADP_add); % uM
options = odeset('NonNegative',[1:p.NOde]); % concentrations of the state variables should be positive

for i=p.ISub:1:p.NSub 
 
X0=ICs(p); p_tem=p; 
%%% solving ODEs and calculating state variables 
T0=-1; %jj=1; %p.Es=1; % extra states 
cc=1;
for ii=1:1:p.ADPL % ii=1: add substrate, ii=2-6 add 4 doses of ADP
    %%% Substrate addition 
    if ii==2
%         X0(p_tem.iMALe)=cc*MAL_index(i)*5e-3; % mM
        X0(p_tem.iMALe)=MAL_index; % M
        X0(p_tem.iPYRe)=cc*PYR_index(i); % M
      
    end

    %%% Solving ODEs
    tspan=[T0:p.tstep:(T0+p.time2(ii,1))];
    [T,X] = ode15s(@ODEs, tspan, X0,options, p_tem);     
    T0=T(end,:);       X0=X(end,:);  % redefining initial values for the next time periode 
    Tc(ii,i)={T};      Xc(ii,i)={X};  CO2(ii,i)=1e3*X(end,p.iO2m);
    %%% Calculating fluxes
    for zz=1:1:length(T) 
        J(zz,:,i)=Fluxes(X(zz,:),p_tem);        
    end
    Jc(ii,i)={J(1:zz,:,i)};
end % ii time protocol for-loop

% storing cell variables in vectors 
Tv(:,i)=        [Tc{2,i}];
Xv(:,:,i)=  [Xc{2,i}]; % M
Jv(:,:,i)=1e9*  [Jc{2,i}]; % nmol/min/mg mito 
end % i substrate for-loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot settings
Position1= [.5,.5, 9, 3];
text_size= 14;
linewidth= 2;
MarkerSizeErr= 2;
%% Simulation OCR dynamic 
set(figure(1),'Units','inches','Position',Position1,'PaperPosition',Position1);
cl =['r','g','b','m','c'];
p.sp=1;
f = [1 5 5 1.5 1];
for i =p.ISub:1:p.NSub
     subplot(1,3,1)
     plot(Time_Data(1:5,:),PYR_Data(1:5,:),'r*',Time_Data(6:10,:),PYR_Data(6:10,:),'go',Time_Data(11:15,:),PYR_Data(11:15,:),'b<',Time_Data(16:20,:),PYR_Data(16:20,:),'md',Time_Data(21:25,:),PYR_Data(21:25,:),'cp',...
         Tv(p.sp:end,i),10^6*f(i)*Xv(p.sp:end,4,i).*Tv(p.sp:end,i).*(p.Ve),cl(i),'linewidth',linewidth)
%      plot(Tv(p.sp:end,i),Jv(p.sp:end,17,i).*Tv(p.sp:end,i)./0.7,cl(i),'linewidth',linewidth)
title('PYR uptake');
     set(gcf,'color','w'); set(gca,'Fontsize',text_size,'linewidth',linewidth);
     hold on
     ylabel('Mass (nmol/mg protein)')
     xlabel('Time (min)')
     text(0.5,190,'D','Fontsize',text_size,'FontWeight','bold')
     xlim([0 10])
     ylim([0 200])
% ax=gca; ax.XLim=[0 inf];
box off;  
end

for j =p.ISub:1:p.NSub
     subplot(1,3,2)
     plot(Time_Data(1:5,:),MAL_Data(1:5,:),'r*',Time_Data(6:10,:),MAL_Data(6:10,:),'go',Time_Data(11:15,:),MAL_Data(11:15,:),'b<',Time_Data(16:20,:),MAL_Data(16:20,:),'md',Time_Data(21:25,:),MAL_Data(21:25,:),'cp',...
         Tv(p.sp:end,j),10^6*j/2*(Xv(p.sp:end,5,j).*Tv(p.sp:end,j)).*p.Ve,cl(j),'linewidth',linewidth)
%      xlabel('Time (sec)');   
title('MAL uptake');
     set(gcf,'color','w'); set(gca,'Fontsize',text_size,'linewidth',linewidth);
     hold on
%      ylabel('Mass (nmol/mg)')
     xlabel('Time (min)')
     text(0.5,190,'E','Fontsize',text_size,'FontWeight','bold')
     xlim([0 10])
     ylim([0 200])
% ax=gca; ax.XLim=[0 inf];
box off;  
end
for z =p.ISub:1:p.NSub
     subplot(1,3,3)
%      plot(Tv(p.sp:end,z),1.5*10^8*Xv(p.sp:end,6,z).*Tv(p.sp:end,z).*p.Ve./0.5,cl(z),'linewidth',linewidth)
     plot(Time_Data(1:5,:),CIT_Data(1:5,:),'r*',Time_Data(6:10,:),CIT_Data(6:10,:),'go',Time_Data(11:15,:),CIT_Data(11:15,:),'b<',Time_Data(16:20,:),CIT_Data(16:20,:),'md',Time_Data(21:25,:),CIT_Data(21:25,:),'cp',...
         Tv(p.sp:end,z),10^8*5*z*(Xv(p.sp:end,6,z)-Xv(5,6,z)).*p.Ve,cl(z),'linewidth',linewidth)
%      xlabel('Time (sec)');   
title('CIT production');
     set(gcf,'color','w'); set(gca,'Fontsize',text_size,'linewidth',linewidth);
     hold on
%      ylabel('Mass (nmol/mg)')
     xlabel('Time (min)')
     text(0.5,190,'F','Fontsize',text_size,'FontWeight','bold')
     ylim([0 200])
     xlim([0 10])
legend('PYR=0mM','PYR=0.5mM','PYR=1mM','PYR=5mM','PYR=10mM','Fontsize',10)
legend box off
box off;  
end

% xlim([0 35]);  xticks([0:7:35])
% ylim([0 80]);  yticks([0:20:80])

