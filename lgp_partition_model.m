%  LGP_PARTITION_MODEL: Fit conditional densities using a partition model
%                       framework.
%
%  syntax:
%  results = lgp_partition_model(y,X,varargin)
%
%  y: nx1 depedent variable to obtain a density estimate of
%  X: nxp design matrix of covariates, each column representing a
%    different variable.  Analysis will be done on the standardized
%    version (done in code below).
%
%  Optional Arguments:
%  niter: number of MCMC iterations. Default = 10,000
%  burn: number of burn in iterations.  Default = 1,000.
%  m: number of grid points for the density estimation.  Default = 400.
%  Mmax: maximum number of partitions. Defaults to 10.
%  n_min: minimum number of observations required in a partition.  Defaults
%    to 50.
%  precision: a tuning parameter for the proposal distribution for the
%    weights, w.  Higher precision decreases the proposal variance and
%    yields higher acceptance rates of moves in the weight vector w.
%    Default: 1.
%  printyes:  (0 or 1) Whether or not to print statistics throughout MCMC.
%    Defaults to 1.
%  w: a px1 vector of initial values of the weights.  
%     Default is equal weight (1/sqrt(p),...,1/sqrt(p)).
%
%  The output of the function is a structure with the following elements:
%  Mpost:  A vector with the number of partitions on each iteration of
%    the MCMC chain.
%  Spost:  A cell object with the indices of the data points which are
%    centers of the tesselation at each iteration of the MCMC chain.
%  acceptancepercent:  the percentage of time a move was accepted in the
%    MCMC chain (including burnin)
%  wacceptancepercent:  the percentage of time a move was accepted when a
%    change in weights was proposed.
%  W: A niter by p matrix with the weights at each iteration of the MCMC
%  llike: the marginal likelihood of the data on each iteratin of the
%    MCMC algorithm.
%
%  Use lgp_graph to get some basic graphs of the output.  
% 
%  See also, LGP_GRAPH.

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

function results = lgp_partition_model(y,X,varargin)
  
  
  % Input Parser
  ip = inputParser;
  ip.FunctionName = 'lgp_partition_model';
  ip.addOptional('burn',1000);
  ip.addOptional('m',400);
  ip.addParameter('Mmax',0);
  ip.addOptional('n_min',50);
  ip.addOptional('niter',10^4);
  ip.addOptional('precision',[]);
  ip.addOptional('printyes',1);
  ip.addOptional('w',[]);
  % ip.addOptional('xstar',[]);
  ip.addRequired('X');
  ip.addRequired('y');
  
  ip.parse(y,X,varargin{:});
  
  m = ip.Results.m;
  Mmax = ip.Results.Mmax;
  niter = ip.Results.niter;
  n_min = ip.Results.n_min;
  burn = ip.Results.burn;
  precision = ip.Results.precision;
  printyes = ip.Results.printyes;
  w = ip.Results.w;

  % Standardize the covariate space;
  [X,~,~] = zscore(X);
  
  % Check for duplicate values.
  %   If duplicate values, add small amount of jitter to those values
  [C,~,IC] = unique(X,'rows');
  if(size(C,1) < size(X,1))
      warning('Duplicate values of X discovered.  Adding random jitter.')
      thetab = tabulate(IC);
      thetab = thetab(thetab(:,2) > 1,:);
      jitterindex = ismember(IC,thetab(:,1));
      newX = X;
      newX(jitterindex,:) = newX(jitterindex,:) + normrnd(0,10^(-9),sum(jitterindex),size(X,2));
      X = newX;
  end
 
  % Set original weights (if not supplied)
  p = size(X,2);
  if size(w,1) < p
      w = ones(p,1)/sqrt(p);
      disp('NOTE: Starting MCMC with equal covariate (tesselation) weights');
  end
  
  if isempty(precision)
     precision = 1;
     disp('NOTE: No precision specified.  Defaulting to 1.');
  end
    
  n = size(y,1);
  % Assign Mmax if not assigned
  if Mmax == 0
      Mmax = 10;
      disp('NOTE: Maximum number of partitions set to 10')
  end
  
  xmin = min([min(y) mean(y) - 3*std(y)]);
  xmax = max([max(y) mean(y) + 3*std(y)]);
  Z = linspace(xmin,xmax,m);

  % Set up the general GP structure
  % Priors for the sigma2 and l (Assuming the grid points have been centered
  %   and scaled)
  pm = prior_sqrtt('s2',10,'nu',4);
  pl = prior_t('s2', 1, 'nu', 4);
  % Default Parameters 
  sigma2 = 1;
  % Best guess of lengthscale (From gpsmooth function in lgpdens from
  %    Riihimaki's excellent MATLAB code)
  Xn = zscore(Z);
  h=max(diff(Xn(end,1:2)).^2,1/length(y).^(1/5)/2);
  
  % With Prior
  % NOTE: the starting value for lengthScale was obtained by looking in the
  %   gpsmooth function in the lgpdens function from Riihimaki.
  cf = gpcf_sexp('magnSigma2',sigma2,'magnSigma2_prior',pm,...
      'lengthScale_prior',pl,'lengthScale',h*repmat(2,[1 size(Xn,1)]));
  gpmflin = gpmf_linear('prior_mean',0,'prior_cov',100);
  % NOTE: I have turned interactions off since this was not mentioned in
  %    the original LGP density paper for 1D (probably doesn't matter for 1D).
  gpmfsq = gpmf_squared('prior_mean',0,'prior_cov',100,'interactions','off');
  % Set GP
  gpgeneral = gp_set('lik',lik_lgp,'cf',cf,'jitterSigma2',1e-6,'meanf',{gpmflin,gpmfsq},...
      'latent_method', 'Laplace');

  %%%%%%%%%%%%%%%%%
  %%%%% MCMC %%%%%%
  %%%%%%%%%%%%%%%%%
  
  % Acceptance Percent for w
  waccept = 0;
  wtotal = 0;
  
  % Initial values  
  S = randsample(n,1); % Index of tesselation centers (random start point)
  M = size(S,1); % number of partitions
  prt = multpartfunc_w(X,S,w); % Current Tesselation Index
  % Log-likelihood
  [llike, GP] = llikefunc(gpgeneral,prt,y,M,Z);
  
  % Initialize vectors which hold posterior samples
  W = zeros(niter,p); % Weights
  Mpost = zeros(niter,1); % Number of partitions
  Spost = cell(niter,1); % Index on X, listing tesselation centers
  LLIKE = zeros(niter,1); % Keep track of marginal log-likelihood on each iteration.
  
  naccepted = 0; % counting acceptance rate of the MCMC chain
  display('Starting MCMC...');
  for ii=1:(niter+burn)
    accepted = 0;
    % Determine which step to do 
    if M > 1 && M < Mmax
        r = randsample(4,1);
        birthrev = 0;
        deathrev = 0;
    elseif M == 1
        r = randsample([1 3 4],1);
    elseif M == Mmax
        r = randsample([2 3 4],1);
    else
        error('A case that was not expected occurred in MCMC')
    end
    
   if r == 1 % Birth
      % Choose an index not already in S
      ind = ismember(1:n,S);
      IND = 1:n;
      IND = IND(ind == 0);
      Sprop = [S; randsample(IND,1)]; % Proposed Tesselation Center
      Mprop = M + 1;
      % proposed partition
      prtprop = multpartfunc_w(X,Sprop,w);
      % Make sure we have a minimum number of obs in each partition
      tab = tabulate(prtprop);
      if(all(tab(:,2) > n_min)) % Proceed only if we have enough obs in each partition
          % Calculate log likelihood
          [llikeprop, GPprop] = llikefunc(gpgeneral,prtprop,y,Mprop,Z);
          % Corrections on boundary cases to the log MH ratio
          if M > 1 && M < (Mmax - 1)
              birthrev = 0;
          elseif M == 1
              birthrev = log(3/4);
          elseif M == (Mmax - 1)
              birthrev = log(4/3);
          else
              error('Unexpected case in Birth Step')
          end
          % log-MH ratio
          lr = llikeprop - llike + birthrev;
          if lr > log(rand(1))
              llike = llikeprop;
              S = Sprop;
              M = Mprop;
              prt = prtprop;
              accepted = 1;
              GP = GPprop;
          end
      end
    
    elseif r == 2 % Death Step
        % Randomly delete index from existing centers
        Sprop = S;
        Sprop(randsample(size(S,1),1)) = [];
        Mprop = M - 1;
        % proposed partition
        prtprop = multpartfunc_w(X,Sprop,w);

        % Calculate log likelihood
        [llikeprop, GPprop] = llikefunc(gpgeneral,prtprop,y,Mprop,Z);
        % Corrections on boundary cases to the log MH ratio
        if M > 2 && M < Mmax
            deathrev = 0;
        elseif M == 2
            deathrev = log(4/3);
        elseif M == Mmax
            deathrev = log(3/4);
        else
            error('Unexpected case in Death Step')
        end
        % log-MH ratio
        lr = llikeprop - llike + deathrev;
        if lr > log(rand(1))
            llike = llikeprop;
            S = Sprop;
            M = Mprop;
            prt = prtprop;
            accepted = 1;
            GP = GPprop;
        end
    
    elseif r == 3 % Move Step
      mind = randsample(M,1);
      Sprop = S;
      Sprop(mind) = []; % Delete one index
      % Add in one index from ones not previously in S
      ind = ismember(1:n,S);
      IND = 1:n;
      IND = IND(ind == 0);
      Sprop = [Sprop; randsample(IND,1)];
      
      % Proposed M
      Mprop = M;
      
      % proposed partition
      prtprop = multpartfunc_w(X,Sprop,w);
      % Make sure we have a minimum number of obs in each partition
      tab = tabulate(prtprop);
      if(all(tab(:,2) > n_min)) % Proceed only if we have enough obs in each partition
          % Calculate log likelihood
          [llikeprop, GPprop] = llikefunc(gpgeneral,prtprop,y,Mprop,Z);
          % No boundary cases to correct for with a move step

          % log-MH ratio
          lr = llikeprop - llike;
          if lr > log(rand(1))
              llike = llikeprop;
              S = Sprop;
              M = Mprop;
              prt = prtprop;
              accepted = 1;
              GP = GPprop;
          end
      end
    elseif r == 4 % Change weights step
        wtotal = wtotal + 1;
               
        % Propose a new weight using a Dirichlet distribution as the
        %   proposal distribution
        % Precision parameter (higher precision increases acceptance rates)
        alpha = precision*w.^2;
        
        wprop = sqrt(drchrnd(alpha',1)');
        alphaprop = precision*wprop.^2;
        
        prtprop = multpartfunc_w(X,S,wprop);
        % Make sure we have a minimum number of obs in each partition
        tab = tabulate(prtprop);
        if(all(tab(:,2) > n_min)) % Proceed only if we have enough obs in each partition
            % Calculate log likelihood
            [llikeprop, GPprop] = llikefunc(gpgeneral,prtprop,y,M,Z);
            % log-MH ratio
            % Uniform prior means we don't have anything but the likelihood 
            %   to determine acceptance since we have uniform priors on
            %   everything
            % Add in the dirichlet penalty since it is not "symmetric"
            pen = drchpdf(w'.^2,alphaprop',1) - drchpdf(wprop'.^2,alpha',1);
            lr = llikeprop - llike + pen;

            if lr > log(rand(1))
                llike = llikeprop;
                w = wprop;
                prt = prtprop;
                accepted = 1;
                waccept = waccept + 1;
                GP = GPprop;
            end
        end
    end
  
    naccepted = naccepted + accepted;
  
    if printyes
      if ii <= burn && mod(ii,10) == 0
        disp(['Burn-in ',num2str(ii),'/',num2str(burn),...
            ', log-lik = ',num2str(llike),...
            ', acceptances = ',num2str(naccepted/ii),...
            ' & ',num2str(waccept/wtotal),...
            ', Partitions = ',num2str(M)]) % print every 10 iterations
      elseif ii > burn && mod(ii,10) == 0
          disp(['Iteration ',num2str(ii-burn),'/',num2str(niter),...
              ', log-lik = ',num2str(llike),...
              ', acceptance = ',num2str(naccepted/ii),...
              ' & ',num2str(waccept/wtotal),...
              ', Partitions = ',num2str(M)]) % print every 10 iterations
      end
    end
    
    if ii > burn
      % Now Store Values
      Mpost(ii-burn) = M;
      Spost{ii-burn} = S;
      W(ii-burn,:) = w;
      LLIKE(ii-burn) = llike;
    end
  end
    
  % Acceptance rates
  acceptancepercent = naccepted/(niter+burn);
  wacceptancepercent = waccept/wtotal;

  % Function Return
  results = struct('Mpost',Mpost,...
      'Spost',{Spost},'acceptancepercent',acceptancepercent,...
      'wacceptancepercent',wacceptancepercent,'W',W,...
      'llike',LLIKE);
    
  % End of Function
  
  
  
  function [ll, GP] = llikefunc(gp,prt,y,M,Z)
  % This function calculates the marginal likelihood, but also returns
  %   the gaussian processes with their optimized hyperparameters in GP.
  %   GP can then be passed to another function (mdleavgdraw) to 
  %   obtain the mean and covariance matrix of the latent function f.
  ll = 0;
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
  
    gp=gp_optim(gp,Zscaled,nypart,'opt',opt, 'optimf', @fminlbfgs);

    % Change GP structure to drop the priors for the 'l' and 'sigma2'
    %   The LGP paper just optimizes the hyperparameters and then
    %   just uses those values.  The priors for the hyperparameters are
    %   simply to aid in finding the smoothness parameters.  The marginal 
    %   distribution in the paper is given the MAP hyperparameter values.
    gp.cf{1} = gpcf_sexp(gp.cf{1},'magnSigma2_prior',[],'lengthScale_prior',[]);
    GP{ii} = gp;
    
    % Calculate Marginal likelihood and add it to the rest.
    ll = ll - gpla_e([],gp,'x',Zscaled,'y',nypart);
  end
  
  
  
  
  
