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

% This file analyzes posterior output and produces figures
%   for supplementary_E.m

addpath(genpath('../src'))
addpath(genpath('../gpstuff'))
% Read in data;
dat = readtable('../data/supplementary_E.csv','Delimiter',',','ReadVariableNames',false);
y = dat{:,1};
X = dat{:,2};

load('../output/supplementary_E_rep1/mcmc_id1.mat');
    
[~,I] = max(output.llike + output.lprior);
plot(output.llike + output.lprior)
Xn = zscore(X);
prt = multpartfunc_w(Xn,output.Spost{I},output.W(I,:)');
figure()
scatter(X(:,1),y,50,prt,'filled');

% 100% correct classification
min(X(prt == 1,1))
max(X(prt == 1,1))
min(X(prt == 2,1))
max(X(prt == 2,1))
min(X(prt == 3,1))
max(X(prt == 3,1))

figure()
subplot(1,3,2)
lgpdens(y(prt == 3), 'percentiles', [2.5, 97.5]);
hold on
    xx = linspace(0,1,1000);
    plot(xx,betapdf(xx,10,30),'--','LineWidth',2)
hold off
xlim([0,1]);
ylim([0,10]);
title('.25 < X < .5');

subplot(1,3,3)
lgpdens(y(prt == 2),'bounded',[1 1],'range',[0 1], 'percentiles', [2.5, 97.5])
hold on
    xx = linspace(.01,.99,1000);
    plot(xx,betapdf(xx,.5,.5),'--','LineWidth',2)
hold off
xlim([0,1]);
ylim([0,10]);
title('X > .5');

subplot(1,3,1)
lgpdens(y(prt == 1), 'percentiles', [2.5, 97.5]);
hold on
    xx = linspace(0,1,1000);
    plot(xx,betapdf(xx,30,20),'--','LineWidth',2)
hold off
xlim([0,1]);
ylim([0,10]);
title('X < .25');

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 7 3];
%saveas(fig,'../figs/supplementary_Epdens.jpg')
%print('../figs/supplementary_Edens','-depsc')

% Save data and make official paper figures in R
[p1, pq1, xt1] = lgpdens(y(prt == 1), 'percentiles', [2.5, 97.5]);
[p2, pq2, xt2] = lgpdens(y(prt == 2), 'percentiles', [2.5, 97.5], 'bounded',[1 1],'range',[0 1]);
[p3, pq3, xt3] = lgpdens(y(prt == 3), 'percentiles', [2.5, 97.5]);
p = [p1, p2, p3];
pq = [pq1, pq2, pq3];
xt = [xt1, xt2, xt3];
rdat = [xt, p, pq];
csvwrite('Rfigs/data/supplementary_E.csv', rdat);