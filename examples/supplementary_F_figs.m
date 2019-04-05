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
%   for supplementary_F_figs.m

addpath(genpath('../src'))
addpath(genpath('../gpstuff'))
% Read in data;
dat = readtable('../data/supplementary_F.csv','Delimiter',',','ReadVariableNames',false);
y = dat{:,1};
X = dat{:,2};

load('../output/supplementary_F_rep1/mcmc_id1.mat');
    
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



lgp_graph(output,y,X)


figure()
subplot(1,3,2)
lgpdens95(y(prt == 1))
hold on
    xx = linspace(0,1,1000);
    plot(xx,betapdf(xx,10,30),'--','LineWidth',2)
hold off
xlim([0,1]);
ylim([0,10]);
title('.25 < X < .5');

subplot(1,3,3)
lgpdens95(y(prt == 2),'bounded',[1 1],'range',[0 1])
hold on
    xx = linspace(.01,.99,1000);
    plot(xx,betapdf(xx,.5,.5),'--','LineWidth',2)
hold off
xlim([0,1]);
ylim([0,10]);
title('X > .5');


subplot(1,3,1)
lgpdens95(y(prt == 3))
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
saveas(fig,'../figs/supplementary_Fpdens.jpg')
print('../figs/supplementary_Fdens','-depsc')




threes = output.Mpost == 3;
[~,I] = max(output.llike(threes) + output.lprior(threes));
threesnum = find(threes);
I2 = threesnum(I);
prt2 = multpartfunc_w(Xn,output.Spost{I2},output.W(I2,:)');
scatter(X(:,1),X(:,2),50,prt2,'filled');
line([0,.75],[0,.75]); line([.75,1],[.75,.75]); line([.75,.75],[.75,1])
output.llike(I2) + output.lprior(I2)

[ll1,allLL1,gp1,gpgen1,lprior1] = troubleshoot(y,X,output.Spost{II},output.W(II,:),output.Mpost(II));
[ll1_wprior,allLL1_wprior,gp1_wprior,gpgen1_wprior,lprior1_wprior] = troubleshoot(y,X,output.Spost{II},output.W(II,:),output.Mpost(II));
[ll2,allLL2,gp2,gpgen2,lprior2] = troubleshoot(y,X,output.Spost{I2},output.W(I2,:),output.Mpost(I2));

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
