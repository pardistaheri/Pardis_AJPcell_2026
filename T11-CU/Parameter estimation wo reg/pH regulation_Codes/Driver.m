
close all; clear; clc
% warning off; format short e

%% Set initial guesses and bounds for the parameters to be optimized
          
K_Hi = 1e-07; %Proton 1st binding constant (M)
Tmax = 1.0411289e-07; %Maximum forward reaction rate(mmol/min)
nh =3;

p0 = [K_Hi, Tmax,nh];  %Binding Constants initial guess

%% Set the optimizer settings and run the optimizer for parameter estimation
% %% Set the optimizer settings and run the optimizer for parameter estimation
options = optimset('fmincon');
options = optimset(options,...
    'TolFun',1e-8,...
    'TolX',1e-8,...
    'Display','iter',...
    'Maxiter',500,...
    'MaxFunEvals',5000);

numIters = 100;
mpars = []; fvals = [];
randLb = 0.5; randUb = 1.5;
numPars = length(p0);

for i = 1:numIters
    temp = p0;
    randNum = (randLb + (randUb-randLb)*rand(1,numPars));
    %Generates random number b/w 0.5 and 1.5 effectively allowing p0 to
    %vary between 0.5p0 and 1.5p0 (a maximal increase or decrease of 50%) 

    p0 = p0.*randNum;
    lbp = p0/1000; ubp = p0*1000;
    [mpar,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@Error,p0,[],[],[],[],lbp,ubp,[],options);
    
    mpars(i,:) = mpar;
    fvals(i) = fval;
    p0 = temp;
end
[minFval,minIdx] = min(fvals);
mpar = mpars(minIdx,:);
minFval
mpar
save('mpar.txt','mpar','-ascii');

disp('Estimated parameter values'); disp(mpar);
disp('Function at optimization termination'); disp(fval);
disp('Exitflag at optimization termination'); disp(exitflag);
disp('Output at optimization termination'); disp(output);
disp('Lambda at optimization termination'); disp(lambda);
disp('Gradient at optimization termination'); disp(grad);
disp('Hessian at optimization termination'); disp(hessian);

%% Run simulations and plot results
 load mpar.txt
 run Simulations;

