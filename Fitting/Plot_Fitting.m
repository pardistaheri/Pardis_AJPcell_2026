%% Sequential ADP figures Normoxia
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% declare variables & read parameters 
% X: state variables;  T: time in min;  J: fuxes;  Pk: peak values of fluxes;  HPk: half the peak values  
% c stands for cell structure; v stand for vector 
run Read_Condt
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax


% Read data
JO2_Pkd=[data.JO2_Pk(:,1) data.JO2_Pk(:,2) data.JO2_Pk(:,3) data.JO2_Pk(:,4)]; % nM
JO2_Pk_Err=[data.JO2_Pk_Err(:,1) data.JO2_Pk_Err(:,2) data.JO2_Pk_Err(:,3) data.JO2_Pk_Err(:,4)]; % nM
JO2_nmol=[data.JO2_nmol_TC(:,1) data.JO2_nmol_TC(:,2) data.JO2_nmol_TC(:,3) data.JO2_nmol_TC(:,4)]; % nM
dPsi=[data.dPsi(:,1) data.dPsi(:,2) data.dPsi(:,3)]; % nM
RCI_data=[data.RCI(:,1) data.RCI(:,2) data.RCI(:,3) data.RCI(:,4)]; % nM
RCI_Err=[data.RCI_Err(:,1) data.RCI_Err(:,2) data.RCI_Err(:,3) data.RCI_Err(:,4)]; % nM

    % substrates and ADP additions 
    PYR_index = [1 0 0 0];
    GLU_index = [0 1 0 0];
    MAL_index = [1 1 0 0];
    SUC_index = [0 0 1 1];
    ADP_add = [12.5 25 50 100] * 1e-6;
    p.ADPL = length(ADP_add);
    
    options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
    
    
    % Preallocate cell arrays
    Tc = cell(p.ADPL + 3, p.NSub); % p.Es = 3
    Xc = cell(p.ADPL + 3, p.NSub);
    Jc = cell(p.ADPL + 3, p.NSub);
    Pkc = cell(p.ADPL + 3, p.NSub);
    HPkc = cell(p.ADPL + 3, p.NSub);
    RCRc = cell(p.ADPL + 3, p.NSub);
    
    for i = p.ISub:1:p.NSub 
        X0 = ICs(p); 
        p_tem = p; 
        
        T0 = 0; %start of time
        jj = 1; %ADP addition counter
        Es = 3; % extra states (state 1, state 2, and state 4) 
        cc = 1;
        
        % Initialize Pk array for this substrate
        Pk = zeros(p.ADPL + Es, 1);
        HPk = zeros(p.ADPL + Es, 1);
        RCR = zeros(p.ADPL + Es, 1);
        
        for ii = 1:(p.ADPL + Es)
            %%% Substrate addition 
            if i == 4 % Suc+Rot, inhibition of ComplexI from the start of experiment
                p_tem.ini_VTmax(p.iCI) = 0 * p_tem.ini_VTmax(p.iCI);
            end 
            
            if ii == (Es - 1)
                X0(p_tem.iMALe) = cc * MAL_index(i) * 5e-3;
                X0(p_tem.iPYRe) = cc * PYR_index(i) * 10e-3;
                X0(p_tem.iGLUe) = cc * GLU_index(i) * 10e-3;
                X0(p_tem.iSUCe) = cc * SUC_index(i) * 7e-3;
            end
            
            %%% ADP addition
            if ii >= Es && ii ~= (p.ADPL + Es)
                X0(p_tem.iADPe) = X0(p.iADPe) + ADP_add(jj);
                jj = jj + 1;
            end
            
            %%% Solving ODEs
            tspan = T0:p.tstep:(T0 + p.time(ii,1));
            [T, X] = ode15s(@ODEs, tspan, X0, options, p_tem);  
            
            T0 = T(end);        
            X0 = X(end,:);
            
            Tc{ii,i} = T;
            Xc{ii,i} = X;
            
            %%% Calculating fluxes - FIXED HERE
            n_fluxes = length(Fluxes(X(1,:), p_tem)); % Get number of fluxes
            J = zeros(length(T), n_fluxes); % Preallocate J with correct dimensions
            
            for zz = 1:length(T)
                flux_result = Fluxes(X(zz,:), p_tem);
                J(zz,:) = flux_result; % This should work now
            end
            
            Jc{ii,i} = J;
            
            s = 14; % CIV flux- OCR
            st = 20; 
            
            %%% finding peaks of state 3 & RCI
            if ii < Es 
                Pk(ii,i) = mean(J(st:end-20, s));
                HPk(ii,i) = mean(J(st:end-20, s));
                RCR(ii,i) = 1;
            else
                Pk(ii,i) = max(J(:, s));
                HPk(ii,i) = 0.5 * Pk(ii,i);
                if ii > 2 % Make sure we have a valid denominator
                    RCR(ii,i) = Pk(ii,i)*1.35 / Pk(2,i);
                else
                    RCR(ii,i) = 1;
                end
            end
            
            Pkc{ii,i} = Pk(ii,i);
            HPkc{ii,i} = HPk(ii,i);
            RCRc{ii,i} = RCR(ii,i);
        end
        
        % storing cell variables in vectors 
Tv(:,i)=        [Tc{1,i};   Tc{2,i};    Tc{3,i};    Tc{4,i};    Tc{5,i};    Tc{6,i};]; % min
Xv(:,:,i)=1e0*  [Xc{1,i};   Xc{2,i};    Xc{3,i};    Xc{4,i};    Xc{5,i};    Xc{6,i};]; % M
Jv(:,:,i)=1e9*  [Jc{1,i};   Jc{2,i};    Jc{3,i};    Jc{4,i};    Jc{5,i};    Jc{6,i};]./1.4; % nmol/min/mg mito
JPkv(:,:,i)=1e9*[Pkc{1,i}; Pkc{2,i};  Pkc{3,i};  Pkc{4,i};  Pkc{5,i};  Pkc{6,i};]./1.4; % nmol/min/mg mito
RCRv(:,:,i)=[RCRc{1,i}; RCRc{2,i}; RCRc{3,i}; RCRc{4,i}; RCRc{5,i}; RCRc{6,i};]; %unitless
end


%%%%%%%%%%%Figure 1
Position1= [.25,.25, 3, 9];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];

set(figure(1),'Units','inches','Position',Position1,'PaperPosition',Position1)
 

for i=p.ISub:1:p.NSub
    subplot(3,1,1)
plot(Tv(:,i),JO2_nmol(:,i),cl(i),'linestyle','-.','linewidth',2);

set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
hold on
     
ylabel('OCR (nmol/min/mg)');
xlabel('Time (min)');
ylim([0 50])
title('Control','Fontsize',text_size2,'FontWeight','bold')
text(1,50,'A','Fontsize',text_size2,'FontWeight','bold')
box off;
end

for i=p.ISub:1:p.NSub 
% plots
%%% CIV fluxes 
subplot(3,1,2)
p.sp=1; % start of plot 

OCR(:,i)=movmean(Jv(p.sp:end,p.iCIV,i),5); % smooth the flux

plot(Tv(p.sp:end,i),OCR(:,i),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
hold on
legend('P+M','G+M','Suc','Suc+Rot','Fontsize',text_size1,'FontWeight','bold')     
ylabel('OCR (nmol/min/mg)');
xlabel('Time (min)');
text(1,50,'D','Fontsize',text_size2,'FontWeight','bold')
legend box off
ylim([0 50])
box off;
end


for i=p.ISub:1:p.NSub 

subplot(3,1,3)
Pkvx=[0 12.5 25 50 100];
plot(Pkvx,JPkv(2:6,:,i),cl(i),'linewidth',2); hold on

errorbar(Pkvx, JO2_Pkd(:,i), JO2_Pk_Err(:,i),cl(i),'Marker','o','MarkerSize',MarkerSizeErr,...
        'MarkerFaceColor',cl(i),'LineStyle','none','linewidth',linewidth);
    
scatter(Pkvx,JO2_Pkd(:,i),cl(i),'linewidth',1);

set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
     hold on

% ylim([0 700]*1e0); 
ylabel('OCR (nmol/min/mg)');
xlabel('ADP conc. (\muM)');
ylim([0 50])
text(5,50,'G','Fontsize',text_size2,'FontWeight','bold')
box off;
end

% for i=p.ISub:1:p.NSub
%     Pkvx=[0 12.5 25 50 100];
% subplot(4,1,4)
% plot(Pkvx,RCRv(2:6,:,i),cl(i),'linewidth',2); hold on
% 
%     
% errorbar(Pkvx, RCI_data(:,i), RCI_Err(:,i),cl(i),'Marker','o','MarkerSize',MarkerSizeErr,...
%         'MarkerFaceColor',cl(i),'LineStyle','none','linewidth',linewidth);
%     
% scatter(Pkvx,RCI_data(:,i),cl(i),'linewidth',1);
% 
% set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
% ylabel('RCR','Fontsize',12);
% xlabel('ADP conc. (\muM)');
% ylim([0 8])
% text(5,8,'J','Fontsize',text_size2,'FontWeight','bold')
% box off;
% end

%%%%%%%%%%%Figure 4
Position4= 3*[.5,.5, 5.5, 5];
text_size1= 10;
text_size2= 12;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;

set(figure(4),'Units','inches','Position',Position4,'PaperPosition',Position4)
titles = ["PDH","CITS","ACON","IDH","AKGDH","SCAS","NDK","FH",...
    "MDH","GOT","CI","CII","CIII","CIV","CV","AK","PYRH","GLUH","DCCS",...
    "DCCM","TCC","OME","GAE","ANT","PIC","HLEAK"];
p.sp=1;
for i =p.ISub:1:p.NSub
   for j= 1:1:26
     subplot(7,4,j)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,j,i),cl(i),'linewidth',linewidth)
%      xlabel('Time (sec)');   

title(titles(j));
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     hold on
     xlim([3 20])
     ylabel('J(nmol/min/mg)')
% ax=gca; ax.XLim=[0 inf];
box off;  
   end
end

%%%%%%%%%%% Figure 5
set(figure(5),'Units','inches','Position',Position4,'PaperPosition',Position4);


titles = ["ADPe","ATPe","Pie","PYRe","MALe","CITe","AKGe","SUCe",...
    "GLUe","ASPe","He","ADPm","ATPm","AMPm","GDPm","GTPm","NADm","NADHm","UQm",...
    "UQH2m","CytCox","CytCred","Pim","GLUm","ASPm","PYRm",...
    "OXAm","CITm","ICITm","AKGm","SCOAm","SUCm","FUMm",...
    "MALm","COAm","ACOAm","O2m","dPsim","Hm"];


for i =p.ISub:1:p.NSub
   for j= 1:1:39
     subplot(7,6,j)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,j,i),cl(i),'linewidth',linewidth)
%      xlabel('Time (sec)');   
if j==7 
    ylim([0 10^-11]);
end
if j==38 
    ylim([140 180]);
end
title(titles(j));
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     hold on
     ylabel('Conc.(M)')
     xlim([2 20])
% ax=gca; ax.XLim=[0 inf];
box off; 
   
end
end

%%%%%%%%%%%Figure 6
Position6= [.5,.5, 3, 3];
text_size1= 10;
text_size2= 12;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];

set(figure(6),'Units','inches','Position',Position6,'PaperPosition',Position6)


for i=p.ISub:1:p.NSub
scatter(Tv(:,i),JO2_nmol(:,i),cl(i)); hold on
OCR(:,i)=movmean(Jv(p.sp:end,p.iCIV,i),5); % smooth the flux
plot(Tv(p.sp:end,i),OCR(:,i),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
hold on
     
ylabel('OCR (nmol/min/mg)');
xlabel('Time (min)');
ylim([0 50])
% legend('P+M','G+M','Suc','Suc+Rot','Fontsize',text_size2,'FontWeight','bold')
title('Normoxia','Fontsize',14,'FontWeight','bold')
text(1,50,'A','Fontsize',text_size2,'FontWeight','bold')
% legend box off
box off;
end

%%%%%%%%%%% Figure 7
Position7= [.5,.5, 3, 15];
text_size1= 10;
text_size2= 12;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
set(figure(7),'Units','inches','Position',Position7,'PaperPosition',Position7);
 
for i =3:1:p.NSub
    subplot(5,1,1)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,32,i).*10^6,cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('SUC_m(\muM)','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     title('Normoxia','Fontsize',14,'FontWeight','bold')
     legend('Suc','Suc+Rot','Fontsize',text_size2,'FontWeight','bold')
     ylim([0 30000])
     legend box off
     box off
     subplot(5,1,2)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,33,i).*10^6,cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('MAL_m(\muM)','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     ylim([0 500])
     box off;
     subplot(5,1,3)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,27,i).*10^6,cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('OXA_m(\muM)','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     ylim([0 0.001])
     box off;
     subplot(5,1,4)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,23,i).*10^6,cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('Pi_m(\muM)','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     ylim([1 1.5]*10^4)
     box off;
     subplot(5,1,5)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,38,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('\Delta\Psi(mV)','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')  
     ylim([100 200])
     box off;
end
 
%%%%%%%%%%% Figure 8
Position8= [.5,.5, 3, 15];
text_size1= 10;
text_size2= 12;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
set(figure(8),'Units','inches','Position',Position8,'PaperPosition',Position8);
 
for i =3:1:p.NSub
    subplot(5,1,1)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,19,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_D_C_C_S','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     title('Normoxia','Fontsize',14,'FontWeight','bold')
     legend('Suc','Suc+Rot','Fontsize',text_size2,'FontWeight','bold')
     ylim([-20 60])
     legend box off
     box off
     subplot(5,1,2)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,20,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_D_C_C_M','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     ylim([-50 0])
     box off;
     subplot(5,1,3)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,12,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_S_D_H','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     ylim([0 50])
     box off;
     subplot(5,1,4)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,25,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_P_I_C','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     ylim([-10 70])
     box off;
     subplot(5,1,5)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,14,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_O_2','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')  
     ylim([0 50])
     box off;
end
  
