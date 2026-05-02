function X = ICs(p)

% Set initial conditions (concentration are in molar unless otherwise noted)
c=1e0;
%% Buffer space 
% -- Nucleotides --
IC(p.iADPe) = 0;
IC(p.iATPe) = 0;

% -- Substrates --
IC(p.iPie) = c*5e-3;
IC(p.iPYRe) = 0;
IC(p.iMALe) = 0;
IC(p.iCITe) = 0;
IC(p.iaKGe) = 0;
IC(p.iSUCe) = 0;    
IC(p.iGLUe) = 0; 
IC(p.iASPe) = 0; % 10

% -- Ions --
IC(p.iHe) = 10^(-7.15);

%% Mito matrix space 
% -- redox pools --
Total_TCA_Pool = c*1e-3; %0.1e-3;
Total_NADH = c*3e-3;
Total_CytC = c*3e-3;
Total_ANP = c*10e-3; 
Total_COA = c*1e-3; % 0.925e-3
Total_UQ = c*1.5e-3;

% -- Nucleotides --
IC(p.iADPm) = 0.025*Total_ANP;
IC(p.iATPm) = 0.975*Total_ANP;
IC(p.iAMPm) = 0;
IC(p.iGDPm) = 0.1*IC(p.iADPm); % 15
IC(p.iGTPm) = 0.1*IC(p.iATPm); 

% -- respiratory substrates --
IC(p.iNADm) = 0.25*Total_NADH;
IC(p.iNADHm) = 0.75*Total_NADH;
IC(p.iUQm) = 0.25*Total_UQ;
IC(p.iUQH2m) = 0.75*Total_UQ;

IC(p.iCytCoxi) = 0.5*Total_CytC; % 20
IC(p.iCytCred) = 0.5*Total_CytC;
IC(p.iPim) = 1*10e-3; % 22 % 10 by Dr Dash

% -- TCA substrates --
IC(p.iGLUm) = Total_TCA_Pool/10;
IC(p.iASPm) = Total_TCA_Pool/10; % 25
IC(p.iPYRm) = Total_TCA_Pool/10;   
IC(p.iOXAm) = Total_TCA_Pool/10;
IC(p.iCITm) = Total_TCA_Pool/10;
%IC(p.iICITm) = Total_TCA_Pool/10;
IC(p.iICITm) = 0;
IC(p.iaKGm) = Total_TCA_Pool/10; % 30
IC(p.iSCOAm) = Total_TCA_Pool/10;
IC(p.iSUCm) = Total_TCA_Pool/10; 
IC(p.iFUMm) = Total_TCA_Pool/10;
IC(p.iMALm) = Total_TCA_Pool/10;
IC(p.iCOAm) = 0.5*(Total_COA); % 35
IC(p.iACOAm) = 0.5*(Total_COA);
IC(p.iO2m) = 210e-6; 
IC(p.idPsi) = 140; %mV % 38
%% -- ions --
IC(p.iHm) = 10^(-7.6); % 39
IC(p.iMgm)= 0.8e-3;
IC(p.iCam) = 100e-9; 
%% All initial conditions
X = IC;
