%     This file is part of bayes-cde.
% 
%     bayes-cde is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     bayes-cde is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with bayes-cde.  If not, see <http://www.gnu.org/licenses/>.
%
%     Copyright 2017-2019 Richard Payne.

% This file provides posterior analysis and plots for the 
%  piecewise constant model in section4_2.m

addpath(genpath('../src'))
addpath(genpath('../gpstuff'))
% Read in data;
dat = readtable('../data/section4_2.csv','Delimiter',',','ReadVariableNames',false);
y = dat{:,1};
X = dat{:,2:3};

if 0 
    % Get MAP estimate from posterior samples
    load('../output/section4_2_rep1/mcmc_id1.mat');
    output1 = output;
    load('../output/section4_2_rep2/mcmc_id1.mat');
    output2 = output;
    load('../output/section4_2_rep3/mcmc_id1.mat');
    output3 = output;
    [lpost1,I1] = max(output1.llike + output1.lprior);
    [lpost2,I2] = max(output2.llike + output2.lprior);
    [lpost3,I3] = max(output3.llike + output3.lprior);
    [~,III] = max([lpost1,lpost2,lpost3]);
    output = eval(strcat(['output',num2str(III)]));
    [~,I] = max(output.llike + output.lprior);
    Xn = zscore(X);
    prt = multpartfunc_w(Xn,output.Spost{I},output.W(I,:)');
    % plot data and posterior partition
    scatter(X(prt==1,1),X(prt==1,2),50,prt(prt==1),'ko');
    hold on
        scatter(X(prt==2,1),X(prt==2,2),50,prt(prt==2),'r+');
        scatter(X(prt==3,1),X(prt==3,2),50,prt(prt==3),'bv');
    hold off
    line([0,.75],[0,.75],'Color','black'); 
    line([.75,1],[.75,.75],'Color','black'); 
    line([.75,.75],[.75,1],'Color','black');
    title('Posterior Partition');
    xlabel('x_1');
    ylabel('x_2');
    % saveas(gcf,'../figs/section4_2.jpg')
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 4 4];
    print('../figs/section4_2','-depsc')
     
    % plot densities in each partition
    lgp_graph(output,y,X)
    
    % plot traceplot of log-posterior
    plot(output.llike + output.lprior)
    
    % Plot posterior vs true densities   
    figure()
    subplot(1,3,1);
    lgpdens(y(prt == 1), 'percentiles', [2.5 97.5])
    hold on
        xx = linspace(0,11,500);
        plot(xx,gampdf(xx,10,1/2),'--','LineWidth',2)
    hold off
    title('X_1 < X_2, X_2 < .75')
    ylabel('density')
    subplot(1,3,2);
    lgpdens(y(prt == 2), 'percentiles', [2.5 97.5])
    hold on
        xx = linspace(-3,9,500);
        plot(xx,.5*normpdf(xx,1,1) + .5*normpdf(xx,5,1),'--','LineWidth',2)
    hold off
    xlabel('y');
    title('X_1 > X_2, X_1 < .75')
    subplot(1,3,3);
    lgpdens(y(prt == 3), 'percentiles', [2.5 97.5])
    hold on
        xx = linspace(-2,4,500);
        plot(xx,normpdf(xx,1,sqrt(.5)),'--','LineWidth',2)
    hold off
    title('X_1 > .75, X_2 > .75')
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    % saveas(fig,'../figs/section4_2dens.jpg')
    print('../figs/section4_2dens','-depsc') 
    
    % MCMC Diagnostics
    plot([output1.llike + output1.lprior;
         output2.llike + output2.lprior;
         output3.llike + output3.lprior])
    title('Log Posterior')
    xlabel('Iteration')
    % number of partitions during MCMC
    plot([output1.Mpost; output2.Mpost; output3.Mpost])
    % traceplots of weights
    plot([output1.W; output2.W; output3.W])
        
    % Look at the tempered chain
    load('../output/section4_2_rep1/mcmc_id8.mat');
    output18 = output;
    load('../output/section4_2_rep2/mcmc_id8.mat');
    output28 = output;
    load('../output/section4_2_rep3/mcmc_id8.mat');
    output38 = output;
    % acceptance percents for birth, death, change, change weights
    (output18.all_acc_percs*20000 + ...
        output28.all_acc_percs*40000 + ...
        output38.all_acc_percs*40000) ./ 100000
        
    % Look at the tempered chains
    SIZES = [];
    swapaccepts = [];
    LPOST = zeros(10^5,8);
    biglpostmax = [];
    biglpostmin = [];
    TOTS = zeros(8,4);
    ACC = TOTS;
    for ii=1:8
        load(strcat('../output/section4_2_rep1/mcmc_id',num2str(ii),'.mat'));
        output1tmp = output;
        load(strcat('../output/section4_2_rep2/mcmc_id',num2str(ii),'.mat'));
        output2tmp = output;
        load(strcat('../output/section4_2_rep3/mcmc_id',num2str(ii),'.mat'));
        output3tmp = output;
        LPOST(:,ii) = [output1tmp.llike + output1tmp.lprior;
            output2tmp.llike + output2tmp.lprior;
            output3tmp.llike + output3tmp.lprior];
%         ACCEPTS(ii,:) = (output1tmp.all_acc_percs*20000 + ...
%             output2tmp.all_acc_percs*40000 + ...
%             output3tmp.all_acc_percs*40000) ./ 100000;
        TOTS(ii,:) = output1tmp.n_proposed + ...
            output2tmp.n_proposed + ...
            output3tmp.n_proposed;
        ACC(ii,:) = output1tmp.n_accepted + ...
            output2tmp.n_accepted + ...
            output3tmp.n_accepted;
        %SIZES(ii).counts = tabulate([output1tmp.Mpost; output2tmp.Mpost; output3tmp.Mpost]);
        tmpcounts = tabulate([output1tmp.Mpost; output2tmp.Mpost; output3tmp.Mpost]);
        SIZES(ii,1:size(tmpcounts,1)) = tmpcounts(:,3);
        % Find largest log posterior for each partition size
        for jj=1:size(tmpcounts,1)
            Mpost = [output1tmp.Mpost; output2tmp.Mpost; output3tmp.Mpost];
            ind = Mpost == jj;
            tmplpostmax = max(LPOST(ind,ii));
            tmplpostmin = min(LPOST(ind,ii));
            biglpostmax(ii,jj) = tmplpostmax;
            biglpostmin(ii,jj) = tmplpostmin;
        end
        
        swapaccepts(ii) = (20000*output1tmp.swap_accept + ...
            40000*output2tmp.swap_accept + 40000*output3tmp.swap_accept)/100000;
    end
    biglpostmax(biglpostmax == 0) = NaN;
    biglpostmin(biglpostmin == 0) = NaN;
    ACCEPTS = ACC ./ TOTS;
    sum(sum(ACC)) / 10^5 % Effective search rate
      
    tabulate([output1tmp.Mpost; output2tmp.Mpost; output3tmp.Mpost])
    xnames = {'1','2','3','4'};
    ynames = {'1.0','.98','.94','.85','.67','.43','.22','.10'};
    heatmap(xnames,ynames,SIZES / 10^5)
    xlabel('Partition Size');
    ylabel('Inverse Temperature');
    title('Proportion of Partition Size by Inverse Temperature');
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 4 3];
    saveas(fig,'../figs/section4_2partsize.jpg')
    
    plot(str2double(ynames),swapaccepts*100,'k')
    title('Swap Acceptance %')
    xlabel('Inverse Temperature')
    ylabel('Swap Acceptance %')
    xlim([0,1.1]);
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 4 3];
    saveas(fig,'../figs/section4_2swaps.jpg')
    
    plot(LPOST(:,1),'k')
    xlabel('Iteration');
    ylabel('Log Posterior');
    title('Log-Posterior Traceplot')
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    saveas(fig,'../figs/section4_2trace.jpg')
    print('../figs/section4_2trace','-depsc')
    
    
    heatmap(biglpostmax)
    figure()
    heatmap(biglpostmin)
    
    % Plot distributions of different partition sizes
    for ii=1:8
        load(strcat('../output/section4_2_rep1/mcmc_id',num2str(ii),'.mat'));
        output1tmp = output;
        load(strcat('../output/section4_2_rep2/mcmc_id',num2str(ii),'.mat'));
        output2tmp = output;
        load(strcat('../output/section4_2_rep3/mcmc_id',num2str(ii),'.mat'));
        output3tmp = output;
        figure()
        for jj=1:size(biglpostmax,2)
            Mpost = [output1tmp.Mpost; output2tmp.Mpost; output3tmp.Mpost];
            ind = Mpost == jj;
            sdat = LPOST(ind,ii);
            if jj == 2
                hold on;
            end
            if ~isempty(sdat)
                ksdensity(sdat);
            end
            title(strcat(['ID = ',num2str(ii)]));
        end
        hold off;            
    end
    
    
    
    
    
    load('../output/section4_2_rep1/mcmc_id8.mat');
    output18 = output;
    load('../output/section4_2_rep2/mcmc_id8.mat');
    output28 = output;
    load('../output/section4_2_rep3/mcmc_id8.mat');
    output38 = output;
    (output18.all_acc_percs*20000 + ...
        output28.all_acc_percs*40000 + ...
        output38.all_acc_percs*40000) ./ 100000
    
    load('../output/section4_2_rep1/mcmc_id5.mat');
    output15 = output;
    load('../output/section4_2_rep2/mcmc_id5.mat');
    output25 = output;
    load('../output/section4_2_rep3/mcmc_id5.mat');
    output35 = output;
    (output15.all_acc_percs*20000 + ...
        output25.all_acc_percs*40000 + ...
        output35.all_acc_percs*40000) ./ 100000
    
    
    
    
    
    
    
    threes = output.Mpost == 3;
    [~,I] = max(output.llike(threes) + output.lprior(threes));
    threesnum = find(threes);
    I2 = threesnum(I);
    prt2 = multpartfunc_w(Xn,output.Spost{I2},output.W(I2,:)');
    scatter(X(:,1),X(:,2),50,prt2,'filled');
    line([0,.75],[0,.75]); line([.75,1],[.75,.75]); line([.75,.75],[.75,1])
    output.llike(I2) + output.lprior(I2)
    
    [ll1,allLL1,gp1,gpgen1] = troubleshoot(y,X,output.Spost{II},output.W(II,:),output.Mpost(II));
    [ll2,allLL2,gp2,gpgen2] = troubleshoot(y,X,output.Spost{I2},output.W(I2,:),output.Mpost(I2));
    
    X(output.Spost{II},:)
    X(output.Spost{I2},:)
    
    allLL1(2) + allLL1(4)
    allLL2(2)
    
    % Here's the issue -- it's with the other datasets...
    % Optimization issue?
    allLL1(1) 
    allLL2(1)
    
    % Same here!
    allLL1(3)
    allLL2(3)
    
    allLL1(1) + allLL1(3)
    allLL2(1) + allLL2(3)
    
    
    
    y1_1 = y(prt == 1);
    y1_2 = y(prt2 == 1);
    
    y3_1 = y(prt == 3);
    y3_2 = y(prt2 == 3);
    
    lgpdens(y1_1)
    figure()
    lgpdens(y1_2)
    
    
    y2 = y(prt == 2);
    y4 = y(prt == 4);
    
    mean(y2)
    mean(y4)
    range(y2)
    range(y4)
    histogram(y2,'Normalization','probability'); 
    hold on;  
        histogram(y4,'Normalization','probability'); 
        histogram([y2; y4],'Normalization','probability'); 
        xx = linspace(.01,10,100);
        plot(xx,gampdf(xx,10,1/2));
    hold off;
    
    lgpdens(y2); hold on; lgpdens(y4); hold off;
    
    

    result = lgp_partition_model_newprior(y,X,'niter',10000,'burn',10000,'seed',1,...
        'filepath','../output/sim1_parallel_newprior/');
    load('../output/sim1_parallel_newprior/mcmc_id1.mat')
    lgp_graph2(output,y,X)
    plot(output.W(:,1))
    plot(output.W(:,2))
    output.W(length(output.llike),:)

    result = lgp_partition_model_newprior(y,X,'niter',10000,'burn',10000,'seed',1,...
        'filepath','../output/trashme/','wprop',.01,'uprop',.05);
end
