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
    scatter(X(prt==2,1),X(prt==2,2),50,prt(prt==2),'k+');
    scatter(X(prt==3,1),X(prt==3,2),50,prt(prt==3),'kv');
hold off
line([0,.75],[0,.75],'Color','black'); 
line([.75,1],[.75,.75],'Color','black'); 
line([.75,.75],[.75,1],'Color','black');
% title('Posterior Partition');
xlabel('X_1');
ylabel('X_2');
% saveas(gcf,'../figs/section4_2.jpg')
%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 3.5 3.5];
%print('figs/section4_2','-depsc')

% Plot posterior vs true densities  
figure()
subplot(1,3,1);
lgpdens_dashed(y(prt == 1))
hold on
    xx = linspace(0,11,500);
    plot(xx,gampdf(xx,10,1/2),'k','LineWidth',2)
hold off
title('X_1 < X_2, X_2 < .75')
ylabel('Density')

subplot(1,3,2);
lgpdens_dashed(y(prt == 2))
hold on
    xx = linspace(-3,9,500);
    plot(xx,.5*normpdf(xx,1,1) + .5*normpdf(xx,5,1),'k','LineWidth',2)
hold off
xlabel('Y');
title('X_1 > X_2, X_1 < .75')

subplot(1,3,3);
lgpdens_dashed(y(prt == 3))
hold on
    xx = linspace(-2,4,500);
    plot(xx,normpdf(xx,1,sqrt(.5)),'k','LineWidth',2)
hold off
title('X_1 > .75, X_2 > .75');
%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 10 3.5];
%print('figs/section4_2dens','-depsc') 

% save data for use in other programs (paper graphics made in R due to
%   Matlab graphics issues)
[p1, pq1, xt1] = lgpdens(y(prt == 1), 'percentiles', [2.5, 97.5]);
[p2, pq2, xt2] = lgpdens(y(prt == 2), 'percentiles', [2.5, 97.5]);
[p3, pq3, xt3] = lgpdens(y(prt == 3), 'percentiles', [2.5, 97.5]);
p = [p1, p2, p3];
pq = [pq1, pq2, pq3];
xt = [xt1, xt2, xt3];
rdat = [xt, p, pq];
csvwrite('Rfigs/data/section4_2.csv', rdat);

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

tabulate([output1tmp.Mpost; output2tmp.Mpost; output3tmp.Mpost]);
xnames = {'1','2','3','4'};
ynames = {'1.0','.98','.94','.85','.67','.43','.22','.10'};
heatmap(xnames,ynames,SIZES);
colormap gray
xlabel('Partition size');
ylabel('Inverse temperature');
%title('Proportion of Partition Size by Inverse Temperature');
%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 4 3];
%print('figs/section4_2partsize', '-depsc');

plot(LPOST(:,1),'k')
xlabel('Iteration');
ylabel('Log posterior');
%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 6 3];
%print('figs/section4_2trace','-depsc')