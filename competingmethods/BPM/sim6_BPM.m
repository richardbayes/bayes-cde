addpath('/home/grad/richard/Documents/mallick/density3/competingmethods/BPM/');
dat = table2array(readtable('../../data/sim6.csv','Delimiter',',','ReadVariableNames',false));
dat = dat(:,1:3);
%rng(1531)
%ntest = 10;
%datpred = dat(randsample(size(dat,1),ntest),:);
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
if 0 
    load('sim6_BPM.mat')


    % Look at partition with maximum probability
    [~,I] = max(chainstats.LL_store);
    centers{I}

    % Plot likelihood
    plot(chainstats.LL_store)
    % Plot number of tessellations
    plot(chainstats.k_store)

    % Get posteriors with credible intervals
    ngrid = 100;
    xgrid = linspace(-3,10,ngrid);
    nmcmc = size(preds,1);
    postdens = [];
    for jj=1:size(datpred,1) % For every prediction spot
        tmp = zeros(nmcmc,ngrid);
        for ii=1:nmcmc
            tmp(ii,:) = normpdf(xgrid,themeans(ii,jj),sqrt(sigma2(ii)));
        end
        qtls = quantile(tmp,[.025,.975])';
        pmean = mean(tmp);
        postdens(jj).lb = qtls(:,1);
        postdens(jj).ub = qtls(:,2);
        postdens(jj).m = pmean';
    end
    % Now plot 
    figure()
    figlocs = [3 4 1 2];
    for jj=1:size(datpred,1)
        subplot(2,2,figlocs(jj));
        plot(xgrid,postdens(jj).m,'--k','LineWidth',2)
        hold on
            plot(xgrid,postdens(jj).lb,'--k')
            plot(xgrid,postdens(jj).ub,'--k')
        hold off
        title(strcat('X_1=',num2str(datpred(jj,2)),...
                     ', X_2=',num2str(datpred(jj,3))));
        xlabel('Y');
        ylabel('density');
        % Plot true function...
        hold on
            if datpred(jj,2) > datpred(jj,3) && datpred(jj,3) < .75
                plot(xgrid,gampdf(xgrid,10,.5),'r')
            elseif datpred(jj,2) < datpred(jj,3) && datpred(jj,2) < .75
                bifun = @(x) .5*normpdf(x,1,1) + .5*normpdf(x,5,1);
                plot(xgrid,bifun(xgrid),'r')
            else
                plot(xgrid,normpdf(xgrid,1,sqrt(.5)),'r');
            end
        hold off
    end
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    saveas(fig,'../../figs/sim6BPM.jpg')
    print('../../figs/sim6BPM','-depsc')

    plot(chainstats.LL_store)
    plot(chainstats.k_store)
end