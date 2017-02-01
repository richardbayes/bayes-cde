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
%     Copyright 2017 Richard Payne.
%
%     NOTE: Be sure the GPstuff toolbox is included in the Matlab path.
%       GPstuff must be downloaded separately from:
%       http://research.cs.aalto.fi/pml/software/gpstuff/

% This file runs several examples of conditional density
%   estimation using a partition model framework.

% Generate design matrix
rng(110038);
n = 1000;
X = [rand(n,1)*5, normrnd(0,1,n,1), normrnd(0,5,n,1)];

% MCMC diagnostics
niter = 1000;
burn = 1000;

% Case 1 - Regression.
rng(20382)
y1 = X(:,1) + normrnd(0,.5,n,1);
plot(X(:,1),y1,'o')

rng(111);
start = cputime;
result1 = lgp_partition_model(y1,X(:,1),'niter',niter,...
    'burn',burn,'precision',5);
finish = cputime;
seconds1 = finish - start;

lgp_graph(result1,y1,X(:,1))

% Case 2 - Two variables
rng(0293)
y2 = 3*X(:,1) + normrnd(0,10*abs(X(:,2)),n,1);
plot(X(:,1),y2,'o')
plot(X(:,2),y2,'o')

rng(222);
start = cputime;
result2 = lgp_partition_model(y2,X(:,1:2),'niter',niter,...
    'burn',burn,'precision',5);
finish = cputime;
seconds2 = finish - start;

lgp_graph(result2,y2,X(:,1:2))


% Case 3 - Three Variables
rng(8255)
y3 = 3*X(:,1) - 5*X(:,2) + normrnd(0,3*abs(X(:,3)),n,1);
plot(X(:,1),y3,'o')
plot(X(:,2),y3,'o')
plot(X(:,3),y3,'o')

rng(333);
start = cputime;
result3 = lgp_partition_model(y3,X,'niter',niter,...
    'burn',burn,'precision',5);
finish = cputime;
seconds3 = finish - start;

lgp_graph(result3,y3,X)


% Case 4 - Nuisance Variable (Case 2 with all 3 variables)
rng(333);
start = cputime;
result4 = lgp_partition_model(y2,X,'niter',niter,...
    'burn',burn,'precision',5);
finish = cputime;
seconds4 = finish - start;

lgp_graph(result4,y2,X)

% Final Weights
result4.W(niter,:)

