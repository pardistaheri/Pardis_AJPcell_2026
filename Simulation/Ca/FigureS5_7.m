%% Sequential ADP Ca simulation Normoxia
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% declare variables & read parameters 
% X: state variables;  T: time in min;  J: fuxes;  Pk: peak values of fluxes;  HPk: half the peak values  
% c stands for cell structure; v stand for vector 
addpath('/Users/pardis/Library/CloudStorage/OneDrive-mcw.edu/Computational_Model_Pardis/Lung_mito_model_norm_HT_HS/Normoxia/Simulation/Ca/Enzyme');
load Model_Params
p.org=1; % organ 
p.Algt=1; % algorithm 
p.ISub=1; %3 initial substrate 
p.NSub=3; % number of substrates  and last substrate
p.NPar=27; % number of parameters
p.NOde=43; % number of ODEs
p.close_system=1; % 1 simulates a closed system like Oroboros, and 0 simulates an open system like PTI
p.Q10_corr=p.Q10_con.^((p.Tem-p.Tem_Stnd)/10); % Xiao p.Tem and p.Tem_Stnd are temps at which reaction rate are measured and calculated, respectively 
p.tstep = 1/30;   % time step min
p.tstepplot = 1/30;   % time step min
p.beta=.35;%.35;

Data = load('Data.txt'); 
Time_data=Data(:,1);
PM_data=Data(:,2);
GM_data=Data(:,3);
SUC_data=Data(:,4);

p.Ve=0.75e-3; % buffer_Vol/mito_Vol experimet:.2mg/mL
p.time=1*[2, 2, 3, 2, 3, 2, 3, 2, 5, 2, 1 ]'; % Seq ADP time
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax

    % substrates and ADP additions 
    PYR_index = [1 0 0 0];
    GLU_index = [0 1 0 0];
    MAL_index = [1 1 0 0];
    SUC_index = [0 0 1 1];
    ADP_add =[50 0 50 0 50 0 50] * 1e-6;
    Ca_add = [0 0.1 0 0.2 0 0.6 0] * 1e-6;
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
                X0(p_tem.iCae) = X0(p.iCae) + Ca_add(jj);
                jj = jj + 1;
            end
            
            if ii==10 
                X0(p_tem.idPsi) = 0;
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
                J(zz,:) = flux_result; 
            end
            
            Jc{ii,i} = J;
            
            s = 14; % CIV flux- OCR
            st = 20; 
            
            %%% finding peaks of state 3 & RCI
            if ii < Es 
                Pk(ii,i) = mean(J(st:end-20, s));
                HPk(ii,i) = mean(J(st:end-20, s));
                RCR(ii,i) = 1;
                XPk(ii,:,i)=mean(Xc{ii,i});
                Xmin(ii,:,i)=mean(Xc{ii,i});
            else
                Pk(ii,i) = max(J(:, s));
                HPk(ii,i) = 0.5 * Pk(ii,i);
                XPk(ii,:,i)=max(Xc{ii,i});
                Xmin(ii,:,i)=min(Xc{ii,i});
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
Tv(:,i)=        [Tc{1,i};   Tc{2,i};    Tc{3,i};    Tc{4,i};    Tc{5,i};    Tc{6,i}; Tc{7,i}; Tc{8,i}; Tc{9,i}; Tc{10,i};]; % min
Xv(:,:,i)=1e0*  [Xc{1,i};   Xc{2,i};    Xc{3,i};    Xc{4,i};    Xc{5,i};    Xc{6,i}; Xc{7,i}; Xc{8,i}; Xc{9,i};Xc{10,i};]; % M
% Tv(:,i)=        [Tc{1,i};   Tc{2,i};    Tc{3,i};]; % min
% Xv(:,:,i)=1e0*  [Xc{1,i};   Xc{2,i};    Xc{3,i};]; % M
% Jv(:,:,i)=1e9*  [Jc{1,i};   Jc{2,i};    Jc{3,i};    Jc{4,i};    Jc{5,i};    Jc{6,i}; Jc{7,i}; Jc{8,i}; Jc{9,i};Jc{10,i};]./1.4; % nmol/min/mg mito
% JPkv(:,:,i)=1e9*[Pkc{1,i}; Pkc{2,i};  Pkc{3,i};  Pkc{4,i};  Pkc{5,i};  Pkc{6,i}; Pkc{7,i}; Pkc{8,i}; Pkc{9,i};Pkc{10,i};]./1.4; % nmol/min/mg mito
% RCRv(:,:,i)=[RCRc{1,i}; RCRc{2,i}; RCRc{3,i}; RCRc{4,i}; RCRc{5,i}; RCRc{6,i}; RCRc{7,i}; RCRc{8,i}; RCRc{9,i};RCRc{10,i};]; %unitless
    end
%%%%%%%%%%%Figure 3
Position3= [.25,.25, 10, 4];
text_size1= 18;
text_size2= 22;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
set(figure(3),'Units','inches','Position',Position3,'PaperPosition',Position3)
cl =['r','g','b'];

p.sp=1; % start of plot

subplot(1,2,1)
plot(Time_data,PM_data,'r','linewidth',2); hold on
plot(Time_data,GM_data,'g','linewidth',2); hold on
plot(Time_data,SUC_data,'b','linewidth',2); hold on

set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); 
% legend('P+M','G+M','Suc','Fontsize',text_size1,'FontWeight','bold')
text(0.2,170,'A','Fontsize',text_size2,'FontWeight','bold')
ylabel('Experimental Data \Delta\Psi_m (mV)','Fontsize',18,'FontWeight','bold')
xlabel('Time (min)','Fontsize',18,'FontWeight','bold');
ylim([100 170])
xlim([0 25])
% yticks([80 120 160 200])
% legend box off
box off;
% end

% for i=p.ISub:1:p.NSub 

subplot(1,2,2)
% plot(Tv(30:60,1),100,'r','linewidth',2); hold on
% plot(Tv(30:60,2),100,'g','linewidth',2); hold on
% plot(Tv(30:60,3),100,'b','linewidth',2); hold on
plot(Tv(31:end-30,1),Xv(61:end,p.idPsi,1)-8,'r','linewidth',2); hold on
plot(Tv(31:end-30,2),Xv(61:end,p.idPsi,2)-8,'g','linewidth',2); hold on
plot(Tv(31:end-30,3),Xv(61:end,p.idPsi,3)-8,'b','linewidth',2); hold on

set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); 
legend('P+M','G+M','Suc','Fontsize',text_size1,'FontWeight','bold')
% title('Model Simulation','Fontsize',22,'FontWeight','bold')
text(0.2,170,'B','Fontsize',text_size2,'FontWeight','bold')
ylabel('Model Prediction \Delta\Psi_m (mV)','Fontsize',18,'FontWeight','bold')
xlabel('Time (min)','Fontsize',18,'FontWeight','bold');
ylim([100 170])
xlim([0 25])
% yticks([80 120 160 200])
legend box off
box off;
% end

% %%%%%%%%%%%Figure 4
% Position4= 3*[.25,.25, 5.5, 5];
% text_size1= 10;
% text_size2= 18;
% linewidth= 1;
% markersize=2;
% MarkerSizeErr= 1;
% 
% set(figure(4),'Units','inches','Position',Position4,'PaperPosition',Position4)
% titles = ["PDH","CITS","ACON","IDH","AKGDH","SCAS","NDK","FH",...
%     "MDH","GOT","CI","CII","CIII","CIV","CV","AK","PYRH","GLUH","DCCS",...
%     "DCCM","TCC","OME","GAE","ANT","PIC","HLEAK","CU"];
% p.sp=1;
% for i =p.ISub:1:p.NSub
%    for j= 1:1:27
%      subplot(7,4,j)
%      plot(Tv(p.sp:end,i),Jv(p.sp:end,j,i),cl(i),'linewidth',linewidth)
% %      xlabel('Time (sec)');   
% 
% title(titles(j));
%      set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
%      hold on
%      xlim([2.5 25])
%      ylabel('J(nmol/min/mg)','Fontsize',10,'FontWeight','bold')
%      xlabel('Time (min)','Fontsize',10,'FontWeight','bold');
% % ax=gca; ax.XLim=[0 inf];
% box off;  
%    end
% end
% 
% 
% set(figure(5),'Units','inches','Position',Position4,'PaperPosition',Position4)
% titles = ["ADPe","ATPe","Pie","PYRe","MALe","CITe",'AKGe',"SUCe",...
%     "GLUe","ASPe","He","ADPm","ATPm","AMPm","GDPm","GTPm","NADm","NADHm","UQm",...
%     "UQH2m","CytCox","CytCred","Pim","GLUm","ASPm","PYRm",...
%     "OXAm","CITm","ICITm","AKGm","SCOAm","SUCm","FUMm",...
%     "MALm","COAm","ACOAm","O2m","dPsim","Hm","Cam","Cae","Mgm","Mge"];
% p.sp=1;
% for i =p.ISub:1:p.NSub
%    for j= 1:1:41
%      subplot(7,6,j)
%      plot(Tv(p.sp:end,i),Xv(p.sp:end,j,i),cl(i),'linewidth',linewidth)
% %      xlabel('Time (sec)');   
% 
% title(titles(j));
%      set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
%      hold on
%      xlim([2.5 25])
%      ylabel('Conc(M)','Fontsize',10,'FontWeight','bold')
%      xlabel('Time (min)','Fontsize',10,'FontWeight','bold');
% % ax=gca; ax.XLim=[0 inf];
% box off;  
%    end
% end