% α-Ketoglutarate Dehydrogenase parameter estimation without regulation
% AKG + COA + NAD = SCOA + NADH + CO2 
% A - AKG; B - COA; ; C - NAD; D - SCOA; E - NADH; F - CO2;
%
% NOTE: The following code was adopted from the model developed
% by Qi et al 2011 published in BMC Biochemistry
% Data used for parameter estimations were obtained from Smith et al 1974
% and McMinn et al 1977
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

clear all;
close all;
clc;

% Set initial guesses and bounds for the parameters to be optimized
              
         
K_A =8.4441040e-05; %AKG binding constant(M)
K_B = 2.2126868e-06; %COA binding constant(M)
K_C = 1.5012596e-05; %NAD binding constant(M)
K_D = 7.1794262e-06; %SCOA binding constant(M)
K_E =  2.9672777e-07; %NADH binding constant(M)
Vmaxf = 6.2196047e-01; %Maximum forward reaction rate(mmol/min)



p0 = [K_A,K_B, K_C,K_D,K_E, Vmaxf];  %Binding Constants initial guess
%%% Set the optimizer settings and run the optimizer for parameter estimation
options = optimset('fmincon');
options = optimset(options,...
    'TolFun',1e-8,...
    'TolX',1e-8,...
    'Display','iter',...
    'Maxiter',500,...
    'MaxFunEvals',5000);

numIters = 50;
mpars = []; fvals = [];
randLb = 0.5; randUb = 1.5;
numPars = length(p0);

for i = 1:numIters
    temp = p0;
    randNum = (randLb + (randUb-randLb)*rand(1,numPars));
    %Generates random number b/w 0.5 and 1.5 effectively allowing p0 to
    %vary between 0.5p0 and 1.5p0 (a maximal increase or decrease of 50%) 

    p0 = p0.*randNum;
    lbp = p0/100; ubp = p0*100;
    [mpar,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@Error,p0,[],[],[],[],lbp,ubp,[],options);
    
    mpars(i,:) = mpar;
    fvals(i) = fval;
    p0 = temp;
end

[minFval,minIdx] = min(fvals);
mpar = mpars(minIdx,:);
minFval;


save('mpar.txt','mpar','-ascii');

% disp('Estimated parameter values'); disp(mpar);
% disp('Function at optimization termination'); disp(fval);
% disp('Exitflag at optimization termination'); disp(exitflag);
% disp('Output at optimization termination'); disp(output);
% disp('Lambda at optimization termination'); disp(lambda);
% disp('Gradient at optimization termination'); disp(grad);
% disp('Hessian at optimization termination'); disp(hessian);

%% Run simulations and plot results
 load mpar.txt
 run Simulations

