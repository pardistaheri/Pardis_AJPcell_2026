function ETfluxes=Fluxes(x,p)
Vmaxp=p.ini_VTmax.*p.Q10_corr.*p.pest(1:27); % p.Q10_corr is 1.624
RT=p.R_con*p.Tem;
CO2=1.2582e-02;   

%% state variables 
%%% Buffer region 
% -- Nucleotides --
ADPe = x(p.iADPe);
ATPe = x(p.iATPe);
% -- Substrates --
Pie  = x(p.iPie); 
PYRe = x(p.iPYRe);
MALe = x(p.iMALe); % 5
CITe = x(p.iCITe);
AKGe = x(p.iaKGe);
SUCe = x(p.iSUCe);
GLUe = x(p.iGLUe);
ASPe = x(p.iASPe); % 10
% -- Ions --
He=x(p.iHe);    pH_e=-log10(He);

%%% Matrix region 
% -- Nucleotides --
ADPm   = x(p.iADPm);
ATPm = x(p.iATPm);
AMPm = x(p.iAMPm);
GDPm = x(p.iGDPm); % 15
GTPm = x(p.iGTPm); 
% -- respiratory substrates --
NADm  = x(p.iNADm);
NADHm = x(p.iNADHm);
UQm   = x(p.iUQm);
UQH2m =x(p.iUQH2m); % 20
CytCoxi = x(p.iCytCoxi); 
CytCred = x(p.iCytCred);
Pim = x(p.iPim);
% -- TCA substrates --
GLUm = x(p.iGLUm);
ASPm = x(p.iASPm); % 25
PYRm = x(p.iPYRm); 
OXAm = x(p.iOXAm);
CITm = x(p.iCITm);
ICITm = x(p.iICITm);
AKGm = x(p.iaKGm); % 30
SCOAm = x(p.iSCOAm); 
SUCm = x(p.iSUCm); 
FUMm = x(p.iFUMm);
MALm = x(p.iMALm);
COAm= x(p.iCOAm); % 35
ACOAm= x(p.iACOAm); 
O2=x(p.iO2m); % 37
if (O2 <= 1e-6)
    O2 = 1e-6;
end
dPsi=x(p.idPsi); % 38
% -- ions --
Hm=x(p.iHm);    pH_m=-log10(Hm); % 39

Cam = x(p.iCam); 
Cae = x(p.iCae); 
Mgm = x(p.iMgm);
Mge = x(p.iMge);
Nam = x(p.iNam); 
Nae = x(p.iNae); 
Km = x(p.iKm);
Ke = x(p.iKe);
% Mgm = 0.8e-3;Cam = 100e-9; 
% Km = 140e-3; Nam = 10e-3;
dGH = (p.F_con*dPsi+ p.R_con*p.Tem*log(He/Hm));

%% enzymes fluxes- 14
%1----------------------------------------------
J_PDH =Vmaxp(p.iPDH)*PDH(PYRm,COAm,NADm,ACOAm,CO2,NADHm,Mgm,Cam,ATPm,ADPm,pH_m,p); 
%2----------------------------------------------
J_CITS =Vmaxp(p.iCITS)*CITS(ACOAm,OXAm,COAm,CITm,ATPm,ADPm,AMPm,SCOAm,pH_m,p);
%3----------------------------------------------
J_ACON =Vmaxp(p.iACON)*ACON(CITm,ICITm,pH_m,p);
%4----------------------------------------------
J_IDH =Vmaxp(p.iIDH)*IDH(ICITm,NADm,AKGm,NADHm,CO2,ATPm,ADPm,pH_m,p);
%5----------------------------------------------
J_AKGDH =Vmaxp(p.iAKGDH)*AKGDH(AKGm,COAm,NADm,SCOAm,NADHm,CO2,Cam,ATPm,ADPm,pH_m,p);
%6----------------------------------------------
J_SCAS= Vmaxp(p.iSCAS)*SCAS(SCOAm,GDPm,Pim,SUCm,GTPm,COAm,pH_m,p); 
%7----------------------------------------------
J_NDK=Vmaxp(p.iNDK)*NDK(GTPm,ADPm,GDPm,ATPm,pH_m,p);
%8----------------------------------------------
J_FH =Vmaxp(p.iFH)*FH(FUMm,MALm,ATPm,pH_m,p);
%9----------------------------------------------
J_MDH=Vmaxp(p.iMDH)*MDH(MALm,NADm,OXAm,NADHm,ATPm,ADPm,AMPm,pH_m,p);
%10----------------------------------------------
J_GOT =Vmaxp(p.iGOT)*GOT(ASPm,AKGm,GLUm,OXAm,pH_m,p);
%11----------------------------------------------
J_CI =Vmaxp(p.iCI)*CI(NADHm,UQm,NADm,UQH2m,Cam,dGH,pH_m,p); 
%12----------------------------------------------
J_CII=Vmaxp(p.iCII)*CIISDH(SUCm,UQm,FUMm,UQH2m,OXAm,pH_m,p);
%13----------------------------------------------
J_CIII=Vmaxp(p.iCIII)*CIII(CytCoxi,UQH2m,CytCred,UQm,dPsi,dGH,pH_m,p);
%14------------------------------------------------
J_CIV =Vmaxp(p.iCIV)*CIV(CytCred,O2,CytCoxi,dPsi,dGH,pH_m,p);
%15----------------------------------------------
J_CV = Vmaxp(p.iCV)*CV(ADPm,Pim,ATPm,dGH,pH_m,p);
%16----------------------------------------------
J_AK = Vmaxp(p.iAK)*AK(ADPm,AMPm,ATPm,pH_m,p);

a=1; b=1;
Vmax=b.*[a*Vmaxp(p.iPDH),Vmaxp(p.iCITS),Vmaxp(p.iACON),Vmaxp(p.iIDH),a*Vmaxp(p.iAKGDH),Vmaxp(p.iSCAS),a*Vmaxp(p.iNDK),Vmaxp(p.iFH),...
    a*Vmaxp(p.iMDH),a*Vmaxp(p.iGOT),Vmaxp(p.iCI),Vmaxp(p.iCII),Vmaxp(p.iCIII),a*Vmaxp(p.iCIV),Vmaxp(p.iCV),Vmaxp(p.iAK)];

%% Transporter fluxes- 10
%17----------------------------------------------
T_PYRH=Vmaxp(p.iPYRH)*PYRH(PYRe,He,PYRm,Hm,p);
%18----------------------------------------------------------
T_GLUH=Vmaxp(p.iGLUH)*GLUH(GLUe,He,GLUm,Hm,p);
%19------------------------------------------------
T_DCCS=Vmaxp(p.iDCC1)*DCCS(Pim,SUCe,Pie,SUCm,MALm,p);
%20---------------------------------------------------------------
T_DCCM=Vmaxp(p.iDCC2)*DCCM(Pim,MALe,Pie,MALm,SUCm,CITm,OXAm,Hm,p);
%21-------------------------------------------------------
T_TCC=Vmaxp(p.iTCC)*TCC(CITm,Hm,MALe,CITe,He,MALm,p);
%22-------------------------------------------------------
T_OME=Vmaxp(p.iOME)*OME(AKGe,MALm,AKGm,MALe,p);
%23-----------------------------------------------------------------
T_GAE=Vmaxp(p.iGAE)*GAE(ASPm,GLUe,ASPe,GLUm,Cam,dPsi,p);
%24-------------------------------------------------------------------------
T_ANT =Vmaxp(p.iANT)*ANT(ADPe,ATPm,ADPm,ATPe,dPsi,p);
%25------------------------------------------------------------------------------
T_PIC=Vmaxp(p.iPIC)*PIC(Pie,He,Pim,Hm,p);
%26------------------------------------------------
T_HLEAK=Vmaxp(p.iHLEAK)*HLEAK(He,Hm,dPsi,p);
%27------------------------------------------------
% T_CU=0;
T_CU=Vmaxp(p.iCU)*CU(Cae,Cam,Mge,Mgm,Pie,Pim,dPsi,p);
%28------------------------------------------------
T_NHE=NHE(Nam,He,Nae,Hm,p);
%29------------------------------------------------
T_CHE=CHE(Cam,He,Cae,Hm,p);
%30------------------------------------------------
T_KHE=KHE(Km,He,Ke,Hm,Cam,Mgm,p);
%31------------------------------------------------
T_NCE=NCE(Cam,Nae,Cae,Nam,dPsi,Hm,He,p);

Tmax=b.*[a*Vmaxp(p.iPYRH),Vmaxp(p.iGLUH),Vmaxp(p.iDCC1),Vmaxp(p.iDCC2),Vmaxp(p.iTCC),Vmaxp(p.iOME),...
    Vmaxp(p.iGAE),Vmaxp(p.iANT),a*Vmaxp(p.iPIC),Vmaxp(p.iHLEAK),Vmaxp(p.iCU),];

%% Metabolic reaction and transport fluxes- 24(+TR123)
ETfluxes=[
J_PDH;
J_CITS;
J_ACON;
J_IDH;
J_AKGDH;
J_SCAS;
J_NDK; 
J_FH;
J_MDH;
J_GOT; % 10
J_CI;    
J_CII;
J_CIII;
J_CIV;
J_CV;     % 15
J_AK;
T_PYRH;
T_GLUH;
T_DCCS;
T_DCCM; %20
T_TCC; 
T_OME; 
T_GAE; 
T_ANT;
T_PIC; % 25
T_HLEAK; 
T_CU;
T_NHE;
T_CHE;
T_KHE; % 30
T_NCE;
];