%% clear workspace 
clear all
close all
clc
warning off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read parameters
run Read_Condt
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax
% experimental condition
p.time=[2, 1, 1, 5,1]';
p.NSub=4; % number of substrates  and last substrate
%% substrates and ADP additions 
PYR_index=[1 0 0 0];
GLU_index=[0 1 0 0];
MAL_index=[1 1 0 0];
SUC_index=[0 0 1 1];
ADP_add=[0 100]*1e-6; p.ADPL=length(ADP_add); % uM    
options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
HF =[1:0.1:10]; % Hleak factor
% HF =[0.01 0.02 0.03 0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9 1 1.5]; % Hleak factor
counter1=0; p.Es=3; counter4=20;
counter5=1; % to store T, X and J with different lengths
T0ii=zeros(p.ADPL+p.Es,5);       X0ii=zeros(p.ADPL+p.Es,5);

% Pre-allocate storage for 3D plotting
all_time_data = cell(length(HF), p.NSub);
all_ocr_data = cell(length(HF), p.NSub);
all_dpsi_data = cell(length(HF), p.NSub);
all_o2_data = cell(length(HF), p.NSub);
all_adpm_data = cell(length(HF), p.NSub);
all_atpm_data = cell(length(HF), p.NSub);
all_atpe_data = cell(length(HF), p.NSub);
all_ant_data = cell(length(HF), p.NSub);
all_nadh_data = cell(length(HF), p.NSub);
all_nad_data = cell(length(HF), p.NSub);
all_uqh2_data = cell(length(HF), p.NSub);
all_uq_data = cell(length(HF), p.NSub);
all_cytcred_data = cell(length(HF), p.NSub);
all_cytcoxi_data = cell(length(HF), p.NSub);

for j=1:length(HF)
    HF(j);
for i=p.ISub:1:p.NSub 
X0=ICs(p); p_tem=p; 
% p_tem.pest(p.iCI)=CIF(j)*p.pest(p.iCI); % Proton Leak
p_tem.ini_VTmax(p.iHLEAK)=HF(j)*p_tem.ini_VTmax(p.iHLEAK);
%%% solving ODEs and calculating state variables 
T0=0; jj=1; p.Es=3; % extra states 
cc=1;
for ii=1:1:p.ADPL+p.Es % ii=1: add mito, ii=2: add substrate, ii=3: add Rot i=5, ii=4-9 add 6 doses of ADP
    %%% Substrate addition 
    if i == 4 % Suc+Rot, inhibition of ComplexI from the start of experiment
        p_tem.ini_VTmax(p.iCI) = 0 * p_tem.ini_VTmax(p.iCI);
    end 
    
    
    if ii==p.Es-1
        X0(p_tem.iMALe)=cc*MAL_index(i)*5e-3; % mM
        X0(p_tem.iPYRe)=cc*PYR_index(i)*10e-3; % mM
        X0(p_tem.iGLUe)=cc*GLU_index(i)*10e-3; % mM
        X0(p_tem.iSUCe)=cc*SUC_index(i)*7e-3; % mM
    end
    %%% ADP addition
    if ii>=p.Es && ii~=p.ADPL+p.Es   
        X0(p_tem.iADPe)=X0(p.iADPe)+ADP_add(jj);
        jj=jj+1;
    end
    %%% Solving ODEs
    tspan=[T0:p.tstep:(T0+p.time(ii,1))];
    [T,X] = ode15s(@ODEs, tspan, X0, options, p_tem);      
    T0=T(end,:);       X0=X(end,:);  % redefining initial values for the next time periode 
    TL(ii,i,j)=length(T);
    Tc(ii,i,j)={T};      Xc(ii,i,j)={X}; 
    %%% Calculating fluxes
    for zz=1:1:length(T) %length(tspan) %length(T)
        J(zz,:,i)=Fluxes(X(zz,:),p_tem);        
    end
    Jc(ii,i,j)={J(1:zz,:,i)};

    %%% finding peaks of state 3 for ODEs
     s=14; % CIV flux- OCR
    st=20; 
    %%% finding peaks and length of state 3
   if ii<p.Es 
        Pk(ii,i,j)=mean(J(st:end-20,s,i));  % returns mean of each column   
        HPk(ii,i,j)=mean(J(st:end-20,s,i));    % returns the same value as HPK. We need it for the next for loop
        XPk(ii,i,j)=mean(X(:,i));
        RCR(ii,i,j) = 1;
   elseif ii>=p.Es || ii~=p.ADPL+p.Es
%        if i==3; st=15; else; st=1; end % for state 2 and Suc start from higher values 
        Pk(ii,i,j)=max(J(1:end,s,i)); % returns max of each column 
        HPk(ii,i,j)=.5*max(J(1:end,s,i));
        XPk(ii,i,j)=max(X(:,i));
%         Xmin(ii,i,j)=min(X(:,i));
        RCR(ii,i,j)=Pk(ii,i,j)/Pk(2,i,j);
   end
   
   if ii==4
        Xmin_NADH(j,i)=min(X(:,18));
        Xmean_UQH2(j,i)=mean(X(10:20,20));     
        Xmean_Cytr(j,i)=mean(X(10:20,22));       
        
%         Xmean_NADH(j,i)=mean(X(10:20,18));
%         Xmean_UQH2(j,i)=mean(X(10:20,20));     
%         Xmean_Cytr(j,i)=mean(X(10:20,22));      
        Xmin_dpsi(j,i)= min(X(10:end,38));
        Xmin_ATP(j,i)= min(X(10:end,13));
        Xmax_ATP(j,i)= max(X(:,13));
        
   end

   %%% store variavles in cell format: we use cell because the length of
    %%% arrays might be different for different substrates 
    JPkc(ii,i)={Pk(ii,i)};
    HPkc(ii,i)={HPk(ii,i)};
    RCIc(ii,i)={RCR(ii,i)};
    
end % for ii

% storing cell variables in vectors 
Tv(:,i,j)=      [Tc{2,i,j};    Tc{3,i,j};    Tc{4,i,j};    Tc{5,i,j};];
Xv(:,:,i)=1e0*  [Xc{2,i,j};    Xc{3,i,j};    Xc{4,i,j};    Xc{5,i,j};]; % M
Jv(:,:,i)=1e9*  [Jc{2,i,j};    Jc{3,i,j};    Jc{4,i,j};    Jc{5,i,j};]./1.4;


end
end

Position2= [.25,.25, 3, 15];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
HF_activity =HF*100;
set(figure(14),'Units','inches','Position',Position2,'PaperPosition',Position2)

for i=p.ISub:1:p.NSub 
subplot(5,1,1)
plot(HF_activity,(Xmin_NADH(:,i)./3e-03),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
     
ylabel('NADH ratio','Fontsize',12,'FontWeight','bold');
xlabel('Hleak','Fontsize',12,'FontWeight','bold');
title('Control','Fontsize',18,'FontWeight','bold')
legend('P+M','G+M','Suc','Suc+Rot','Fontsize',text_size1,'FontWeight','bold')
text(120,1,'A','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
xlim([100 900]); xticks([100:200:900])
ylim([0 1]);
legend box off
box off;
end 

for i=p.ISub:1:p.NSub 
subplot(5,1,2)
plot(HF_activity,(Xmean_UQH2(:,i)./1.5e-03),cl(i),'linewidth',2); hold on  
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
ylabel('UQH_2 ratio','Fontsize',12,'FontWeight','bold');
xlabel('Hleak','Fontsize',12,'FontWeight','bold');
xlim([100 900]); xticks([100:200:900])
ylim([0 0.6])
text(120,0.6,'D','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
box off;
end 

for i=p.ISub:1:p.NSub
subplot(5,1,3)
plot(HF_activity,(Xmean_Cytr(:,i)./3e-03),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
ylabel('CytC_r_e_d ratio','Fontsize',12,'FontWeight','bold');
xlabel('Hleak','Fontsize',12,'FontWeight','bold');
xlim([100 900]); xticks([100:200:900])
ylim([0 0.6])
text(120,0.6,'G','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
box off;
end 

for i=p.ISub:1:p.NSub
subplot(5,1,4)
if i==1||i==2
plot(HF_activity,Xmin_dpsi(:,i)-15,cl(i),'linewidth',2); hold on
else 
plot(HF_activity,Xmin_dpsi(:,i)-20,cl(i),'linewidth',2); hold on
end
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on 
ylabel('\Delta\Psi_m (mV)','Fontsize',12,'FontWeight','bold');
xlabel('Hleak','Fontsize',12,'FontWeight','bold');
xtickformat('percentage'); 
xlim([100 900]); xticks([100:200:900])
ylim([100 160])
text(110,160,'J','Fontsize',text_size2,'FontWeight','bold')
box off;
end

for i=p.ISub:1:p.NSub
subplot(5,1,5)
plot(HF_activity,(Xmax_ATP(:,i)-Xmin_ATP(:,i))*10^3,cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on 
ylabel('ATP_m (mM)','Fontsize',12,'FontWeight','bold');
xlabel('Hleak','Fontsize',12,'FontWeight','bold');
xtickformat('percentage'); 
xlim([100 900]); xticks([100:200:900])
ylim([0 4.5])
text(120,4.5,'M','Fontsize',text_size2,'FontWeight','bold')
box off;
end