dat = table2array(readtable('../../data/section4_2.csv','Delimiter',',','ReadVariableNames',false));
dat = dat(:,1:3);
datpred = [.1 .8;
           .8 .1;
           .76 .76;
           .9 .9];
% Y is a dummy variable in this dataset so we can choose the locations
datpred = [ones(size(datpred,1),1),datpred];

% Default Options
STANDARDISE=1; options.standardise=1;
LIN=0; options.linear=LIN;
Near_points=1; options.near_points = Near_points;
k_max = 200; options.k_max = k_max;
mcmc_samples = 80000; options.mcmc_samples = mcmc_samples; 
burn_in = 20000; options.burn_in =burn_in; 
alpha_1=0.1;alpha_2=0.1; options.alpha_1=alpha_1; options.alpha_2=alpha_2; 
SAVE_SAMPLES = 0; options.save=SAVE_SAMPLES; 

rng(55512)
tic
[output,chainstats,preds,themeans,sigma2,centers] = bayes_partition_gauss(dat,datpred,options);
seconds = toc;
save('sim6_BPM.mat')