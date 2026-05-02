%% declare variables 
% X: state variables;  T: time in min;  J: fuxes;  Pk: peak values of fluxes;  HPk: half the peak values  
% c stands for cell structure; v stand for vector 
addpath('/Users/pardis/Library/CloudStorage/OneDrive-mcw.edu/Computational_Model_Pardis/Lung_mito_model_norm_HT_HS/Normoxia/Fitting/Enzyme');
load Model_Params
p.org=1; % organ 
p.Algt=1; % algorithm 
p.ISub=1; % initial substrate 
p.NSub=4; % number of substrates  and last substrate
p.NPar=26; % number of parameters
p.NOde=39; % number of ODEs
p.nH=2.6777;
p.close_system=1; % 1 simulates a closed system like Oroboros, and 0 simulates an open system like PTI
p.Q10_corr=p.Q10_con.^((p.Tem-p.Tem_Stnd)/10); % Xiao p.Tem and p.Tem_Stnd are temps at which reaction rate are measured and calculated, respectively 
p.tstep = 1/30;   % time step min
p.tstepplot = 1/30;   % time step min
p.beta=.35;%.35;
load Exp_data_Lung;  data = data_OM;

p.Ve=2e-3; % buffer_Vol/mito_Vol experimet:.2mg/mL
p.time=1*[2, 3, 2, 3, 4, 6, 2 ]'; % Seq ADP time
p.time1=1*[3, 5, 2]'; % Single ADP time


    