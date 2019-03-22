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

% This file replicates an example in Ma (2017) for use with the
%  logistic Gaussian processes for conditional density estimation.

addpath(genpath('../src'))
addpath(genpath('../gpstuff'))
% Read in data;
dat = readtable('../data/supplementary_F.csv','Delimiter',',','ReadVariableNames',false);
y = dat{:,1};
X = dat{:,2};

% Parallel Example
parpool(8);
result = lgp_partition_model(y,X,'niter',10000,'burn',0,'seed',5315,...
    'filepath','../output/supplementary_F_new_rep1/','saveall',1,'nprint',1000,'n_min',15);

% see supplementary_F_figs.m for posterior analysis


