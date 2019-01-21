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

function [ll, GP, allLL] = llikefunc(gp,prt,y,M,Z)
  % This function calculates the marginal likelihood, but also returns
  %   the gaussian processes with their optimized hyperparameters in GP.
  %   GP can then be passed to another function (mdleavgdraw) to 
  %   obtain the mean and covariance matrix of the latent function f.
  % 
  % Inputs
  % gp: a gp object
  % prt: a vector of integers of the same length of the data where integers
  %      indicate the partition each data point belongs to
  % y: the dependent numeric vector of length n
  % M: The number of partition regions
  % Z: the grid that y is discretized upon for the logistic Gaussian
  %      process.
  %
  % Outputs
  % ll: numeric, the log-likelihood of the partition
  % GP: the optimized gp objects for each partition region
  % allLL: the log-likelihood for each of the M partition regions
  
  ll = 0;
  allLL = zeros(1,M);
  % Calculate and add up the likelihood contribution in each partition
  GP = cell(M,1);
  m = max(size(Z,2));
  gporig = gp;
  % Set optim options
  opt = optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
  for ii=1:M
    % Split y in to the iith partition
    ypart = y(prt == ii);
    
    % Choose a subset of bins (to avoid estimating densities in regions
    %   that don't have any data).
    xmin = min([min(ypart) mean(ypart) - 3*std(ypart)]);
    xmax = max([max(ypart) mean(ypart) + 3*std(ypart)]);
    xlb = find(Z > xmin, 1 ) - 1;
    if xlb < 1
        xlb = 1;
    end
    xup = find(Z < xmax, 1, 'last' ) + 1;
    if xup > m
        xup = m;
    end
    Xpart = Z(xlb:xup);
    
    % Calculate the number in each bin
    nypart = hist(ypart,Xpart);
    nypart = nypart';
    % Optimize Hyperparameters
    
    % Optimization of hyperparameters...
    % Scale X so we can always use the same prior
    % (same strategy is employed in lgpdens just before using gpsmooth
    % function
    Zscaled = ((Xpart - mean(Xpart))./std(Xpart))';
    % Best guess of lengthscale (From gpsmooth function in lgpdens from
    %    riihimaki's excellent MATLAB code)
    h=max(diff(Zscaled(1:2,end)).^2,1/sum(nypart).^(1/5)/2);
    gp = gporig;
    gpcf1 = gp.cf{1};
    gpcf1 = gpcf_sexp(gpcf1, 'lengthScale', h*repmat(2,[1 size(Zscaled,2)]));
    gp = gp_set(gp,'cf',gpcf1);
  
    try 
        gp=gp_optim_rdp(gp,Zscaled,nypart,'opt',opt, 'optimf', @fminlbfgs_rdp);
    catch
        warning('Problem using gp_optim_rdp: returning -Inf for log-likelihood.');
        ll = -Inf;
        return;
    end

    % Change GP structure to drop the priors for the 'l' and 'sigma2'
    %   The LGP paper just optimizes the hyperparameters and then
    %   just uses those values.  The priors for the hyperparameters are
    %   simply to aid in finding the smoothness parameters.  The marginal 
    %   distribution in the paper is given the MAP hyperparameter values.
    gp.cf{1} = gpcf_sexp(gp.cf{1},'magnSigma2_prior',[],'lengthScale_prior',[]);
    GP{ii} = gp;
    
    % Calculate Marginal likelihood and add it to the rest.
    try 
        [~,ll_tmp,~] = gpla_e([],gp,'x',Zscaled,'y',nypart);
        ll_tmp = -ll_tmp;
    catch
        warning('Problem using gpla_e: returning -Inf for log-likelihood.');
        ll = -Inf;
        return;
    end
    allLL(ii) = ll_tmp;
    ll = ll + ll_tmp;
  end
  
  
  
  
  
