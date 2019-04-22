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
% Read in data;
dat = readtable('../data/section4_3.csv','Delimiter',',',...
    'ReadVariableNames',false);
y = dat{:,1};
X = dat{:,2:3};
% load untempered chain (after 20,000 burn in)
load('../output/section4_3_rep2/mcmc_id1.mat');
% find posterior maximum
[~,I] = max(output.llike + output.lprior);
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

% save data for use in other programs (paper graphics made in R due to
%   Matlab graphics issues)
[p1, pq1, xt1] = lgpdens(y(prt == 1), 'percentiles', [2.5, 97.5]);
[p2, pq2, xt2] = lgpdens(y(prt == 2), 'percentiles', [2.5, 97.5]);
[p3, pq3, xt3] = lgpdens(y(prt == 3), 'percentiles', [2.5, 97.5]);
p = [p1, p2, p3];
pq = [pq1, pq2, pq3];
xt = [xt1, xt2, xt3];
rdat = [xt, p, pq];
rdatx = [X, prt];
csvwrite('Rfigs/data/section4_3_densities.csv', rdat);
csvwrite('Rfigs/data/section4_3_X.csv', rdatx);

% Figures similar to paper (R version used for paper)
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
%saveas(fig,'../figs/section4_3.jpg')
%print('../figs/section4_3','-depsc')

plot(1:(8*1e4), output.llike + output.lprior,'k')
ylabel('Log posterior');
xlabel('Iteration');
ylim([-2670, -2643])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 3];
print('figs/section4_3trace','-depsc')