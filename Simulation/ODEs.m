function dxdt = ODEs(~,x,p)
%   The purpose of this function is to establish the rate equations for
%   chemical reactions of molecular compounds and transport of these
%   compounds between compartments in the metabolism model.
ETfluxes=Fluxes(x,p);
%% TCA fluxes
J_PDH=ETfluxes(1);
J_CITS=ETfluxes(2);
J_ACON=ETfluxes(3);
J_IDH=ETfluxes(4);
J_AKGDH=ETfluxes(5);
J_SCAS=ETfluxes(6);
J_NDK=ETfluxes(7);
J_FH=ETfluxes(8);
J_MALD=ETfluxes(9);
J_GOT=ETfluxes(10);
J_CI=ETfluxes(11);
J_CII=ETfluxes(12);
J_CIII=ETfluxes(13);
J_CIV=ETfluxes(14);
J_CV=ETfluxes(15);
J_AK=ETfluxes(16);
%% Transporter fluxes
T_PYRH=ETfluxes(17);
T_GLUH=ETfluxes(18);
T_DCCS=ETfluxes(19);
T_DCCM=ETfluxes(20);
T_TCC=ETfluxes(21);
T_OME=ETfluxes(22);
T_GAE=ETfluxes(23);
T_ANT=ETfluxes(24);
T_PIC=ETfluxes(25);
T_HLEAK=ETfluxes(26);
T_CU=ETfluxes(27);
T_NHE=ETfluxes(28);
T_CHE=ETfluxes(29);
T_KHE=ETfluxes(30);
T_NCE=ETfluxes(31);

%% Net Equations / Mass Balance
%%% Buffer region 
% -- Nucleotides --
dxdt(p.iADPe) = (-T_ANT)/p.Ve;
dxdt(p.iATPe) = (+T_ANT)/p.Ve;
% -- Substrates --
dxdt(p.iPie)= 1*(T_DCCS -T_PIC +T_DCCM)/p.Ve;
dxdt(p.iPYRe) = 1*(-T_PYRH)/p.Ve;
dxdt(p.iMALe) = 1*(-T_DCCM -T_OME -T_TCC)/p.Ve; % 5
dxdt(p.iCITe) = (+T_TCC)/p.Ve;
dxdt(p.iaKGe) = (+T_OME)/p.Ve;
dxdt(p.iSUCe) = (-T_DCCS)/p.Ve;
dxdt(p.iGLUe) = (-T_GLUH -T_GAE)/p.Ve;
dxdt(p.iASPe) = (T_GAE)/p.Ve; % 10
% -- Ions --
dxdt(p.iHe)=0; 
% buffer_m= 2.303*10^(-7.6)/43;  % buffer capacity =0.0043 M in matrix
% dxdt(p.iHe) = -(buffer_m*(T_CHE+T_KHE+T_NHE +T_PYRH +T_GLUH +T_GAE -T_TCC +T_HLEAK +T_PIC ...
%     +J_PDH +2*J_CITS +J_AKGDH +J_SCAS +J_MALD -(4+1)*J_CI -(4-2)*J_CIII -(2+2)*J_CIV +(3-1)*J_CV))/p.Vm; % 39
%%% Mitochondria Martix region 
% -- Nucleotides -
dxdt(p.iADPm) = (T_ANT -J_NDK -J_CV-2*J_AK)/p.Vm;
dxdt(p.iATPm) = (J_CV +J_NDK -T_ANT+J_AK)/p.Vm;
dxdt(p.iAMPm) = (J_AK)/p.Vm;
dxdt(p.iGDPm) = (-J_SCAS +J_NDK)/p.Vm; % 15
dxdt(p.iGTPm) = -dxdt(p.iGDPm);  
% -- respiratory substrates --
dxdt(p.iNADm) = (-J_PDH -J_IDH - J_AKGDH -J_MALD +J_CI)/p.Vm;
dxdt(p.iNADHm) = -dxdt(p.iNADm);
dxdt(p.iUQm) = (-J_CI -J_CII +J_CIII)/p.Vm;
dxdt(p.iUQH2m) = -dxdt(p.iUQm); % 20
dxdt(p.iCytCoxi) = (-2*J_CIII +2*J_CIV)/p.Vi; 
dxdt(p.iCytCred) = -dxdt(p.iCytCoxi);
dxdt(p.iPim)  = (-T_DCCS -T_DCCM -J_SCAS -J_CV +T_PIC)/p.Vm; 
% -- TCA substrates --
dxdt(p.iGLUm) = (J_GOT +T_GAE +T_GLUH)/p.Vm;
dxdt(p.iASPm) = (-J_GOT -T_GAE)/p.Vm; % 25
dxdt(p.iPYRm) = (T_PYRH -J_PDH)/p.Vm; 
dxdt(p.iOXAm) = (-J_CITS +J_MALD+J_GOT)/p.Vm; 
dxdt(p.iCITm) = (-T_TCC +J_CITS -J_ACON)/p.Vm;
dxdt(p.iICITm) = (J_ACON -J_IDH)/p.Vm;
dxdt(p.iaKGm)  = (J_IDH -J_GOT -J_AKGDH -T_OME)/p.Vm; % 30
dxdt(p.iSCOAm)  = (J_AKGDH -J_SCAS)/p.Vm; 
dxdt(p.iSUCm)  = (+T_DCCS +J_SCAS -J_CII)/p.Vm; 
dxdt(p.iFUMm)  =(-J_FH +J_CII)/p.Vm; 
dxdt(p.iMALm)  = (J_FH -J_MALD  +T_DCCM +T_OME +T_TCC)/p.Vm; 
dxdt(p.iCOAm) = (-J_PDH +J_CITS +J_SCAS -J_AKGDH)/p.Vm; % 35
dxdt(p.iACOAm)= (J_PDH-J_CITS)/p.Vm; 
% -- oxygen consumption --
dxdt(p.iO2m) = p.close_system*0.5*(-J_CIV)/(p.Vm +p.Ve); 
%dxdt(p.iO2m) = 250e-06;
% -- memebrane potential --
Cimm=6.75e-6*(p.Vm+p.Vi);   % Mito volume (is 1.3e-6 L in Jason's paper). 
dxdt(p.idPsi) = 1/Cimm*(4*J_CI+ 2*J_CIII+ 4*J_CIV -T_ANT -T_HLEAK -T_GAE -p.nH*J_CV ); % 38
% -- ions / proton leak --
buffer_m= 2.303*10^(-7.6)/43;  % buffer capacity =0.0043 M in matrix
dxdt(p.iHm) = (buffer_m*( T_CHE+ T_KHE +T_NHE+T_PYRH +T_GLUH +T_GAE -T_TCC +T_HLEAK +T_PIC ...
    +J_PDH +2*J_CITS +J_AKGDH +J_SCAS +J_MALD -(4+1)*J_CI -(4-2)*J_CIII -(2+2)*J_CIV +(3-1)*J_CV))/p.Vm; % 39
 
dxdt(p.iCam)=0.0001*(2*T_CU-T_CHE-T_NCE)/p.Vm; %buffering capacity
dxdt(p.iCae)=0.0001*(-2*T_CU+T_CHE+T_NCE)/p.Ve;
dxdt(p.iMgm)=0;
dxdt(p.iMge)=0;
dxdt(p.iNam)=(-T_NHE+3*T_NCE)/p.Vm; 
dxdt(p.iNae)=(T_NHE-3*T_NCE)/p.Ve;
dxdt(p.iKm)=(-T_KHE)/p.Vm; 
dxdt(p.iKe)=(T_KHE)/p.Ve;
%% Output Matrix
dxdt=dxdt';
end

