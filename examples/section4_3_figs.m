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

%     This file analyzes posterior output from the partition model with 
%     data generated from a multivariate normal distribution.

addpath(genpath('../src'))
addpath(genpath('../gpstuff'))
% load untempered chain (after 20,000 burn in)
load('../output/section4_3_rep2/mcmc_id1.mat');
% find posterior maximum
[~,I] = max(output.llike + output.lprior);
% log-posterior traceplot
plot(output.llike + output.lprior)
% normalize data
Xn = zscore(X);
% get partitions for MAP partition
prt = multpartfunc_w(Xn,output.Spost{I},output.W(I,:)');
theinds = output.Spost{I};
% plot partitions by color and the Voronoi centers (black)
figure()
scatter(X(:,1),X(:,2),50,prt,'filled');
hold on
    for ii=1:length(theinds)
        scatter(X(theinds(ii),1),X(theinds(ii),2),50,'k','filled')
    end
hold off
% plot the partition and densities in each partition
lgp_graph(output,y,X)

% Figures for the paper
subplot(2,2,1)
scatter(X(prt==1,1),X(prt==1,2),50,prt(prt==1),'ko');
hold on
    scatter(X(prt==2,1),X(prt==2,2),50,prt(prt==2),'r+');
    scatter(X(prt==3,1),X(prt==3,2),50,prt(prt==3),'bv');
hold off
ylim([3.5,11.5])
xlim([1.5,7.5])
xlabel('x_1');
ylabel('x_2');
Sigma = [1 .5 .1;
         .5 1 .5;
         .1 .5 1];
mu = [1 5 7.5]';
mu1 = mu(1);
mu2 = mu(2:3);
sigma11 = Sigma(1,1);
sigma12 = Sigma(1,2:3);
sigma22 = Sigma(2:3,2:3);    
xx = linspace(-30,30,10000);
for ii=1:max(prt)
    theind = theinds(ii);
    subplot(2,2,ii+1);
    %figure()
    lgpdens(y(prt == ii), 'percentiles', [2.5, 97.5]);
    hold on
        a = X(theind,:)';
        themean = mu1 + sigma12 * (sigma22 \ (a - mu2));
        thecov = sigma11 - sigma12 * (sigma22 \ sigma12');
        plot(xx,normpdf(xx,themean,thecov),'--','LineWidth',2)
        title(strcat(['x_1=',num2str(X(theind,1),2),', x_2=',num2str(X(theind,2),2)]));
    hold off
    xlim([-3,5]);
    ylim([0,.6]);
    ylabel('density');
    xlabel('y');
end
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 6];
saveas(fig,'../figs/section4_3.jpg')
print('../figs/section4_3','-depsc')

figure()
lgpdens(y(prt == 2),'bounded',[1 1],'range',[0 1], 'percentiles', [2.5 97.5]);
hold on
    xx = linspace(0,1,1000);
    plot(xx,betapdf(xx,.5,.5))
hold off

figure()
lgpdens(y(prt == 3), 'percentiles', [2.5 97.5]);
hold on
    xx = linspace(0,1,1000);
    plot(xx,betapdf(xx,30,20))
hold off

% MCMC Diagnostics
% Look at the tempered chains
SIZES = [];
swapaccepts = [];
LPOST = zeros(10^5,8);
biglpostmax = [];
biglpostmin = [];
TOTS = zeros(8,4);
ACC = TOTS;
for ii=1:8
    load(strcat('../output/section4_3_rep1/mcmc_id',num2str(ii),'.mat'));
    output1tmp = output;
    load(strcat('../output/section4_3_rep2/mcmc_id',num2str(ii),'.mat'));
    output2tmp = output;
    LPOST(:,ii) = [output1tmp.llike + output1tmp.lprior;
        output2tmp.llike + output2tmp.lprior];
    TOTS(ii,:) = output1tmp.n_proposed + ...
        output2tmp.n_proposed;
    ACC(ii,:) = output1tmp.n_accepted + ...
        output2tmp.n_accepted;
    tmpcounts = tabulate([output1tmp.Mpost; output2tmp.Mpost]);
    SIZES(ii,1:size(tmpcounts,1)) = tmpcounts(:,3);
    % Find largest log posterior for each partition size
    for jj=1:size(tmpcounts,1)
        Mpost = [output1tmp.Mpost; output2tmp.Mpost];
        ind = Mpost == jj;
        tmplpostmax = max(LPOST(ind,ii));
        tmplpostmin = min(LPOST(ind,ii));
        if isempty(tmplpostmax)
            tmplpostmax = NaN;
        end
        if isempty(tmplpostmin)
            tmplpostmin = NaN;
        end
        biglpostmax(ii,jj) = tmplpostmax;
        biglpostmin(ii,jj) = tmplpostmin;
    end
    swapaccepts(ii) = (20000*output1tmp.swap_accept + ...
        80000*output2tmp.swap_accept)/100000;
end
biglpostmax(biglpostmax == 0) = NaN;
biglpostmin(biglpostmin == 0) = NaN;
ACCEPTS = ACC ./ TOTS;

plot(1:10^5,LPOST(:,1),'k')
title('Log-Posterior Traceplot');
ylabel('Log Posterior');
xlabel('Iteration');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 3];
saveas(fig,'../figs/section4_3_trace.jpg')
print('../figs/sectionr_3_trace','-depsc')