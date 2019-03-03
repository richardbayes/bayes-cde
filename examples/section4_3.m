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

%     This file runs the partition model with data generated from 
%     a multivariate normal distribution.

addpath(genpath('../src'))
addpath(genpath('../gpstuff'))
% Read in data;
dat = readtable('../data/section4_3.csv','Delimiter',',',...
    'ReadVariableNames',false);
y = dat{:,1};
X = dat{:,2:3};

% parallel tempering
parpool(8);
% run MCMC sequentially since saving files is memory intensive
result = lgp_partition_model(y,X,'niter',20000,'burn',0,'seed',7372,...
    'filepath','../output/section4_3_rep1/','saveall',1,'nprint',1000,...
    'n_min',15,'Mmax',20,'sprop_w',.01);
result = lgp_partition_model(y,X,'niter',80000,'burn',0,'seed',37,...
    'filepath','../output/section4_3_rep2/','saveall',1,'nprint',1000,...
    'n_min',15,'Mmax',20,'sprop_w',.01,...
    'resume','../output/section4_3_rep1/');

% See section4_3_figs.m for posterior analysis