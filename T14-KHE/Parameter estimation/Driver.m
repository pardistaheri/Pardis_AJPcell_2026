% KHE binding constants' estimation and biochemical reaction fit
% Km + He -> Ke + Hm
% the flux for KHE transporter in this model was devolped in Audi & Dash lab
% Parameter estimation was done using data from Bazil et al. 2010. the
% The units are as follows Concentration M, Flux mmol/min, Volume =1 ml, Mass in Microg

close all; clear; clc
 warning off; format short e

%% Set initial guesses and bounds for the parameters to be optimized

   
K_A = 1.8e-3; %K Binding constant (M)
K_B = 1.63e-8; %Proton binding constant (M)
Tmax = 1.0896585e-05; %Maximum transport rate(mmol/min)
   
p0 = [K_A, K_B, Tmax];  %Binding Constants initial guess
%%% Set the optimizer settings and run the optimizer for parameter estimation
options = optimset('fmincon');
options = optimset(options,...
    'TolFun',1e-8,...
    'TolX',1e-8,...
    'Display','iter',...
    'Maxiter',500,...
    'MaxFunEvals',5000);

numIters = 20;
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
    fvals(i) = fval
    p0 = temp;
end

[minFval,minIdx] = min(fvals);
mpar = mpars(minIdx,:);
minFval
mpar

save('mpar.txt','mpar','-ascii');

%disp('Estimated parameter values'); disp(mpar);
%disp('Function at optimization termination'); disp(fval);
% disp('Exitflag at optimization termination'); disp(exitflag);
% disp('Output at optimization termination'); disp(output);
% disp('Lambda at optimization termination'); disp(lambda);
% disp('Gradient at optimization termination'); disp(grad);
% disp('Hessian at optimization termination'); disp(hessian);

%% Run simulations and plot results
 load mpar.txt
 run Simulations;

