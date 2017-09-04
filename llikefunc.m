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
  
    gp=gp_optim_rdp(gp,Zscaled,nypart,'opt',opt, 'optimf', @fminlbfgs_rdp);

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
  
  
  
  
  
