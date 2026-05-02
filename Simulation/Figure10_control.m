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
% p.time=[2, 1, 5, 1,1]';
p.time1=1*[3, 5, 2]'; % Single ADP time
p.NSub=3; % number of substrates  and last substrate
% substrates and ADP additions 
PYR_index=[1 0 0 0];
GLU_index=[0 1 0 0];
MAL_index=[1 1 0 0];
SUC_index=[0 0 1 1];
ADP_add=[100]*1e-6; p.ADPL1=length(ADP_add); % uM
options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
CIF =[1:-0.01:0]; % CI factor
counter1=0; p.Es=3; counter4=20;
counter5=1; % to store T, X and J with different lengths
T0ii=zeros(p.ADPL+p.Es,5);       X0ii=zeros(p.ADPL+p.Es,5);

% Pre-allocate arrays for storing results
JPk_CIV = zeros(length(CIF), p.NSub);

for j=1:length(CIF)
%     CIF(j);
    for i=p.ISub:1:p.NSub 
        X0=ICs(p); p_tem=p; 
        p_tem.ini_VTmax(p.iCI)=CIF(j)*p_tem.ini_VTmax(p.iCI);
        
        %%% solving ODEs and calculating state variables 
        T0=-2; jj=1;  
        cc=1;
        
        for ii=1:1:3
            %%% Substrate addition    
    if ii==1
        X0(p_tem.iMALe)=cc*MAL_index(i)*5e-3; % mM
        X0(p_tem.iPYRe)=cc*PYR_index(i)*10e-3; % mM
        X0(p_tem.iGLUe)=cc*GLU_index(i)*10e-3; % mM
        X0(p_tem.iSUCe)=cc*SUC_index(i)*10e-3; % mM
    end
    
    %%% ADP addition
    if ii==2   
        X0(p_tem.iADPe)=X0(p.iADPe)+ADP_add(jj);
        jj=jj+1;     
    end
            
            %%% Solving ODEs
            tspan=[T0:p.tstepplot:(T0+p.time1(ii,1))];
            [T,X] = ode15s(@ODEs, tspan, X0,options, p_tem);     
            T0=T(end,:);       X0=X(end,:);  
            TL(ii,i,j)=length(T);
            Tc(ii,i,j)={T};      Xc(ii,i,j)={X}; 
            
            %%% Calculating fluxes
            for zz=1:1:length(T)
                J(zz,:,i)=Fluxes(X(zz,:),p_tem);        
            end
            Jc(ii,i,j)={J(1:zz,:,i)};
            
            if ii==2
                JPk_CIV(j,i)= 1e9*max(J(:,14,i))./1.4;  % CIV activity
            end
% %             maxJPK_CIV(:,i)=max(JPk_CIV(:,i)); 
        end % for ii
    end % i substrate for-loop
end % j for leak factor


Position2= [.25,.25, 3, 3];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
CI_activity =CIF*100;
set(figure(1),'Units','inches','Position',Position2,'PaperPosition',Position2)


for i=p.ISub:1:p.NSub
plot(CI_activity,JPk_CIV(:,i),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
set(gca, 'XDir', 'reverse');
ylabel('J_O_2 (nmol/min/mg)','Fontsize',12,'FontWeight','bold');
xlabel('CI Activity (%)','Fontsize',12,'FontWeight','bold');
title('Control','Fontsize',18,'FontWeight','bold')
legend('P+M','G+M','Suc','Fontsize',text_size1,'FontWeight','bold')
text(95,45,'A','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
% ytickformat('percentage');
ylim([0 45]); 
box off;
legend box off
end 

%% clear workspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read parameters
run Read_Condt
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax
% experimental condition
p.time1=1*[3, 5, 2]'; % Single ADP time
p.NSub=3; % number of substrates  and last substrate
% substrates and ADP additions 
PYR_index=[1 0 0 0];
GLU_index=[0 1 0 0];
MAL_index=[1 1 0 0];
SUC_index=[0 0 1 1];
ADP_add=[100]*1e-6; p.ADPL1=length(ADP_add); % uM
options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
CIIF =[1:-0.01:0]; % CII factor
counter1=0; p.Es=3; counter4=20;
counter5=1; % to store T, X and J with different lengths
T0ii=zeros(p.ADPL+p.Es,5);       X0ii=zeros(p.ADPL+p.Es,5);

% Pre-allocate arrays for storing results
JPk_CIV = zeros(length(CIIF), p.NSub);

for j=1:length(CIIF)
%     CIF(j);
    for i=p.ISub:1:p.NSub 
        X0=ICs(p); p_tem=p; 
        p_tem.ini_VTmax(p.iCII)=CIIF(j)*p_tem.ini_VTmax(p.iCII);
        
        %%% solving ODEs and calculating state variables 
        T0=-2; jj=1;  
        cc=1;
        
        for ii=1:1:3
            %%% Substrate addition    
    if ii==1
        X0(p_tem.iMALe)=cc*MAL_index(i)*5e-3; % mM
        X0(p_tem.iPYRe)=cc*PYR_index(i)*10e-3; % mM
        X0(p_tem.iGLUe)=cc*GLU_index(i)*10e-3; % mM
        X0(p_tem.iSUCe)=cc*SUC_index(i)*10e-3; % mM
    end
    
    %%% ADP addition
    if ii==2   
        X0(p_tem.iADPe)=X0(p.iADPe)+ADP_add(jj);
        jj=jj+1;     
    end
            
            %%% Solving ODEs
            tspan=[T0:p.tstepplot:(T0+p.time1(ii,1))];
            [T,X] = ode15s(@ODEs, tspan, X0,options, p_tem);     
            T0=T(end,:);       X0=X(end,:);  
            TL(ii,i,j)=length(T);
            Tc(ii,i,j)={T};      Xc(ii,i,j)={X}; 
            
            %%% Calculating fluxes
            for zz=1:1:length(T)
                J(zz,:,i)=Fluxes(X(zz,:),p_tem);        
            end
            Jc(ii,i,j)={J(1:zz,:,i)};
            
            if ii==2
                JPk_CIV(j,i)= 1e9*max(J(:,14,i))./1.4;  % CIV activity
            end
% %             maxJPK_CIV(:,i)=max(JPk_CIV(:,i)); 
        end % for ii
    end % i substrate for-loop
end % j for leak factor


Position2= [.25,.25, 3, 3];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
CII_activity =CIIF*100;
set(figure(2),'Units','inches','Position',Position2,'PaperPosition',Position2)


for i=p.ISub:1:p.NSub
plot(CII_activity,JPk_CIV(:,i),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
set(gca, 'XDir', 'reverse');
ylabel('J_O_2 (nmol/min/mg)','Fontsize',12,'FontWeight','bold');
xlabel('CII Activity (%)','Fontsize',12,'FontWeight','bold');
text(95,45,'D','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage');
ylim([0 45]); 
box off;
end 

%% clear workspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read parameters
run Read_Condt
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax
% experimental condition
p.time1=1*[3, 5, 2]'; % Single ADP time
p.NSub=3; % number of substrates  and last substrate
% substrates and ADP additions 
PYR_index=[1 0 0 0];
GLU_index=[0 1 0 0];
MAL_index=[1 1 0 0];
SUC_index=[0 0 1 1];
ADP_add=[100]*1e-6; p.ADPL1=length(ADP_add); % uM
options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
CIIIF =[1:-0.01:0]; % CI factor
counter1=0; p.Es=3; counter4=20;
counter5=1; % to store T, X and J with different lengths
T0ii=zeros(p.ADPL+p.Es,5);       X0ii=zeros(p.ADPL+p.Es,5);

% Pre-allocate arrays for storing results
JPk_CIV = zeros(length(CIIIF), p.NSub);

for j=1:length(CIIIF)
%     CIF(j);
    for i=p.ISub:1:p.NSub 
        X0=ICs(p); p_tem=p; 
        p_tem.ini_VTmax(p.iCIII)=CIIIF(j)*p_tem.ini_VTmax(p.iCIII);
        
        %%% solving ODEs and calculating state variables 
        T0=-2; jj=1;  
        cc=1;
        
        for ii=1:1:3
            %%% Substrate addition    
    if ii==1
        X0(p_tem.iMALe)=cc*MAL_index(i)*5e-3; % mM
        X0(p_tem.iPYRe)=cc*PYR_index(i)*10e-3; % mM
        X0(p_tem.iGLUe)=cc*GLU_index(i)*10e-3; % mM
        X0(p_tem.iSUCe)=cc*SUC_index(i)*10e-3; % mM
    end
    
    %%% ADP addition
    if ii==2   
        X0(p_tem.iADPe)=X0(p.iADPe)+ADP_add(jj);
        jj=jj+1;     
    end
            
            %%% Solving ODEs
            tspan=[T0:p.tstepplot:(T0+p.time1(ii,1))];
            [T,X] = ode15s(@ODEs, tspan, X0,options, p_tem);     
            T0=T(end,:);       X0=X(end,:);  
            TL(ii,i,j)=length(T);
            Tc(ii,i,j)={T};      Xc(ii,i,j)={X}; 
            
            %%% Calculating fluxes
            for zz=1:1:length(T)
                J(zz,:,i)=Fluxes(X(zz,:),p_tem);        
            end
            Jc(ii,i,j)={J(1:zz,:,i)};
            
            if ii==2
                JPk_CIV(j,i)= 1e9*max(J(:,14,i))./1.4;  % CIV activity
            end
% %             maxJPK_CIV(:,i)=max(JPk_CIV(:,i)); 
        end % for ii
    end % i substrate for-loop
end % j for leak factor


Position2= [.25,.25, 3, 3];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
CIII_activity =CIIIF*100;
set(figure(3),'Units','inches','Position',Position2,'PaperPosition',Position2)


for i=p.ISub:1:p.NSub
plot(CIII_activity,JPk_CIV(:,i),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
set(gca, 'XDir', 'reverse');
ylabel('J_O_2 (nmol/min/mg)','Fontsize',12,'FontWeight','bold');
xlabel('CIII Activity (%)','Fontsize',12,'FontWeight','bold');
text(95,45,'G','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
% ytickformat('percentage'); 
ylim([0 45]); 
box off;
end 

%% clear workspace 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read parameters
run Read_Condt
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax
% experimental condition
p.time1=1*[3, 5, 2]'; % Single ADP time
p.NSub=3; % number of substrates  and last substrate
% substrates and ADP additions 
PYR_index=[1 0 0 0];
GLU_index=[0 1 0 0];
MAL_index=[1 1 0 0];
SUC_index=[0 0 1 1];
ADP_add=[100]*1e-6; p.ADPL1=length(ADP_add); % uM
options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
CVF =[1:-0.01:0]; % CI factor
counter1=0; p.Es=3; counter4=20;
counter5=1; % to store T, X and J with different lengths
T0ii=zeros(p.ADPL+p.Es,5);       X0ii=zeros(p.ADPL+p.Es,5);

% Pre-allocate arrays for storing results
JPk_CIV = zeros(length(CVF), p.NSub);

for j=1:length(CVF)
%     CIF(j);
    for i=p.ISub:1:p.NSub 
        X0=ICs(p); p_tem=p; 
        p_tem.ini_VTmax(p.iCV)=CVF(j)*p_tem.ini_VTmax(p.iCV);
        
        %%% solving ODEs and calculating state variables 
        T0=-2; jj=1;  
        cc=1;
        
        for ii=1:1:3
            %%% Substrate addition    
    if ii==1
        X0(p_tem.iMALe)=cc*MAL_index(i)*5e-3; % mM
        X0(p_tem.iPYRe)=cc*PYR_index(i)*10e-3; % mM
        X0(p_tem.iGLUe)=cc*GLU_index(i)*10e-3; % mM
        X0(p_tem.iSUCe)=cc*SUC_index(i)*10e-3; % mM
    end
    
    %%% ADP addition
    if ii==2   
        X0(p_tem.iADPe)=X0(p.iADPe)+ADP_add(jj);
        jj=jj+1;     
    end
            
            %%% Solving ODEs
            tspan=[T0:p.tstepplot:(T0+p.time1(ii,1))];
            [T,X] = ode15s(@ODEs, tspan, X0,options, p_tem);     
            T0=T(end,:);       X0=X(end,:);  
            TL(ii,i,j)=length(T);
            Tc(ii,i,j)={T};      Xc(ii,i,j)={X}; 
            
            %%% Calculating fluxes
            for zz=1:1:length(T)
                J(zz,:,i)=Fluxes(X(zz,:),p_tem);        
            end
            Jc(ii,i,j)={J(1:zz,:,i)};
            
            if ii==2
                JPk_CIV(j,i)= 1e9*max(J(50:end,14,i))./1.4;  % CIV activity
            end
% %             maxJPK_CIV(:,i)=max(JPk_CIV(:,i)); 
        end % for ii
    end % i substrate for-loop
end % j for leak factor


Position2= [.25,.25, 3, 3];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
CV_activity =CVF*100;
set(figure(4),'Units','inches','Position',Position2,'PaperPosition',Position2)


for i=p.ISub:1:p.NSub
plot(CV_activity,JPk_CIV(:,i),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
set(gca, 'XDir', 'reverse');
ylabel('J_O_2 (nmol/min/mg)','Fontsize',12,'FontWeight','bold');
xlabel('CV Activity (%)','Fontsize',12,'FontWeight','bold');
text(95,45,'J','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
% ytickformat('percentage'); 
ylim([0 45]); 
box off;
end 

%% clear workspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read parameters
run Read_Condt
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax
% experimental condition
p.time1=1*[3, 5, 2]'; % Single ADP time
p.NSub=3; % number of substrates  and last substrate
% substrates and ADP additions 
PYR_index=[1 0 0 0];
GLU_index=[0 1 0 0];
MAL_index=[1 1 0 0];
SUC_index=[0 0 1 1];
ADP_add=[100]*1e-6; p.ADPL1=length(ADP_add); % uM
options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
ANTF =[1:-0.01:0]; % CI factor
counter1=0; p.Es=3; counter4=20;
counter5=1; % to store T, X and J with different lengths
T0ii=zeros(p.ADPL+p.Es,5);       X0ii=zeros(p.ADPL+p.Es,5);

% Pre-allocate arrays for storing results
JPk_CIV = zeros(length(ANTF), p.NSub);

for j=1:length(ANTF)
%     CIF(j);
    for i=p.ISub:1:p.NSub 
        X0=ICs(p); p_tem=p; 
        p_tem.ini_VTmax(p.iANT)=ANTF(j)*p_tem.ini_VTmax(p.iANT);
        
        %%% solving ODEs and calculating state variables 
        T0=-2; jj=1;  
        cc=1;
        
        for ii=1:1:3
            %%% Substrate addition    
    if ii==1
        X0(p_tem.iMALe)=cc*MAL_index(i)*5e-3; % mM
        X0(p_tem.iPYRe)=cc*PYR_index(i)*10e-3; % mM
        X0(p_tem.iGLUe)=cc*GLU_index(i)*10e-3; % mM
        X0(p_tem.iSUCe)=cc*SUC_index(i)*10e-3; % mM
    end
    
    %%% ADP addition
    if ii==2   
        X0(p_tem.iADPe)=X0(p.iADPe)+ADP_add(jj);
        jj=jj+1;     
    end
            
            %%% Solving ODEs
            tspan=[T0:p.tstepplot:(T0+p.time1(ii,1))];
            [T,X] = ode15s(@ODEs, tspan, X0,options, p_tem);     
            T0=T(end,:);       X0=X(end,:);  
            TL(ii,i,j)=length(T);
            Tc(ii,i,j)={T};      Xc(ii,i,j)={X}; 
            
            %%% Calculating fluxes
            for zz=1:1:length(T)
                J(zz,:,i)=Fluxes(X(zz,:),p_tem);        
            end
            Jc(ii,i,j)={J(1:zz,:,i)};
            
            if ii==2
                JPk_CIV(j,i)= 1e9*max(J(:,14,i))./1.4;  % CIV activity
            end
% %             maxJPK_CIV(:,i)=max(JPk_CIV(:,i)); 
        end % for ii
    end % i substrate for-loop
end % j for leak factor


Position2= [.25,.25, 3, 3];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
ANT_activity =ANTF*100;
set(figure(5),'Units','inches','Position',Position2,'PaperPosition',Position2)


for i=p.ISub:1:p.NSub
plot(ANT_activity,JPk_CIV(:,i),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
set(gca, 'XDir', 'reverse');
ylabel('J_O_2 (nmol/min/mg)','Fontsize',12,'FontWeight','bold');
xlabel('ANT Activity (%)','Fontsize',12,'FontWeight','bold');
text(95,45,'M','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
% ytickformat('percentage');
ylim([0 45]); 
box off;
end 