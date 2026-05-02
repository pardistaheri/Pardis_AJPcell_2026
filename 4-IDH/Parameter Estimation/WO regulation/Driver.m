% Isocitrate dehydrogenase parameter estimation without regulation
% ICIT + NAD + H2O = AKG + NADH + CO2 + 2H+
% A - ICIT; B - NAD; ; C - AKG; D - NADH; E - CO2;
%
% NOTE: The following code was adopted from the lung model developed
% by Zhang et al 2018 published in Public Library of Science 
% Data used for parameter estimations were obtained from Qi et al 2008
% published in Biochimica et Biophysica Acta
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

clear all;
close all;
clc;

% Set initial guesses and bounds for the parameters to be optimized

K_A = 4.4715300e-04; %ICIT binding constant(M)
K_B = 4.9382372e-02; %NAD binding constant(M)
K_C = 1.7501750e-04 ; %AKG binding constant(M)
K_D = 6.6718425e-06; %NADH binding constant(M)
Vmaxf = 1.3127534e-02; %Maximum forward reaction rate(mmol/min)



p0 = [K_A,K_B, K_C,K_D, Vmaxf];  %Binding Constants initial guess
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
    lbp = p0/10; ubp = p0*10;
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

