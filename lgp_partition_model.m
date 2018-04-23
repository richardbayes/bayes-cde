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

function output = lgp_partition_model(y,X,varargin)
    % Input Parser
    ip = inputParser;
    ip.FunctionName = 'lgp_partition_model';
    ip.addOptional('burn',1000);
    ip.addOptional('resume',[]);
    addParameter(ip,'filepath','./output/')
    addParameter(ip,'hottemp',.1)
    ip.addOptional('m',400);
    ip.addParameter('Mmax',0);
    ip.addOptional('n_min',50);
    ip.addOptional('niter',10^4);
    ip.addOptional('nprint',10);
    %ip.addOptional('precision',[]);
    ip.addOptional('printyes',1);
    ip.addOptional('saveall',0);
    ip.addOptional('S',[]);
    addParameter(ip,'seed','shuffle');
    ip.addOptional('sprop_w',.1);
    addParameter(ip,'suppress_errors_on_workers',0);
    ip.addOptional('swapfreq',1);
    ip.addOptional('w',[]);
    % ip.addOptional('xstar',[]);
    ip.addRequired('X');
    ip.addRequired('y');

    ip.parse(y,X,varargin{:});

    filepath = ip.Results.filepath;
    hottemp = ip.Results.hottemp;
    m = ip.Results.m;
    Mmax = ip.Results.Mmax;
    niter = ip.Results.niter;
    n_min = ip.Results.n_min;
    nprint = ip.Results.nprint;
    burn = ip.Results.burn;
    resume = ip.Results.resume;
    % precision = ip.Results.precision;
    printyes = ip.Results.printyes;
    S = ip.Results.S;
    saveall = ip.Results.saveall;
    seed = ip.Results.seed;
    sprop_w = ip.Results.sprop_w;
    suppress_errors_on_workers = ip.Results.suppress_errors_on_workers;
    swapfreq = ip.Results.swapfreq;
    w = ip.Results.w;
    
    rng(seed); % Set RNG for the client...
    
    if suppress_errors_on_workers
        disp('NOTE: Errors suppressed on workers.');
    end
    
    % Make output drive
    makefolder_status = mkdir(filepath);
    if ~makefolder_status
        error('Could not make output directory');
    end

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

    
    % Initial values  
    % Set original weights (if not supplied)
    p = size(X,2);
    if isempty(resume)
        if length(w) < p 
            w = ones(p,1)/sqrt(p);
            disp('NOTE: Starting MCMC with equal covariate (tesselation) weights');
        else
            if size(w,1) < size(w,2)
                w = w';
            end
        end

        if isempty(S)
            S = randsample(n,1); % Index of tesselation centers (random start point)
            M = size(S,1); % number of partitions
        else
            if size(S,1) < size(S,2)
                S = S';
            end
            M = size(S,1);
        end


        [prt,pflag] = multpartfunc_w(X,S,w); % Current Tesselation Index
        if pflag
            error('Tessellation has equivalent centers due to the S and w combination. Change start values.')
        end
        % Log-likelihood
        [llike, ~] = llikefunc(gpgeneral,prt,y,M,Z);
        lprior = get_lprior(n,size(X,2),M,Mmax);
    end
    
    % Initialize vectors which hold posterior samples
    W = zeros(niter,p); % Weights
    Mpost = zeros(niter,1); % Number of partitions
    Spost = cell(niter,1); % Index on X, listing tesselation centers
    LLIKE = zeros(niter,1); % Keep track of marginal log-likelihood on each iteration.
    LPRIOR = zeros(niter,1);

    naccepted = 0; % counting acceptance rate of the MCMC chain
    totals = [0 0 0 0];
    accepts = [0 0 0 0];
    swapaccepttotal = 0;
    swaptotal = 0;
    swapaccepttotal_global = 0;
    swaptotal_global = 0;

    % Parallel Tempering Prep
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        poolsize = 0;
    else
        poolsize = poolobj.NumWorkers;
    end
    mm = poolsize;
    % Harmonic Temperatures
    %delta = (1/hottemp - 1)/(mm-1);
    %temps = [1, 1./(1 + delta*(1:mm - 1))];
    % Sigmoidal temperatures
    j1 = log( 1/(-1 + 1.01) - 1);
    jm = log( 1/(-hottemp + 1.01) - 1);
    dm = (jm-j1)/(mm-1);
    js = j1:dm:jm;
    temps = 1.01 - 1./(1 + exp(js));

    spmdsize = min([poolsize,mm]);
    if ~isempty(resume)
        files = dir2(resume);
        Sres = cell(length(files),1);
        wres = cell(length(files),1);
        for ii=1:spmdsize
            try 
                load(strcat([resume,files(ii).name]));
            catch
                error('Cannot load files to resume MCMC. Verify parpool size and number of files match.')
            end
            nnn = length(output.Spost);
            Sres{ii} = output.Spost{nnn};
            wres{ii} = output.W(nnn,:)';
        end
        clear output;
    end
    
        
    if spmdsize < 2
        disp('NOTE: Must have at least two processes to do parallel tempering.');
        disp('      Initiate outside the function with parpool function.');
        % error('Must have at least two processes to do parallel tempering.');
        spmdsize = 1;
    end
    disp(['NOTE: Using ',num2str(spmdsize),' cores.']);
  
    if(spmdsize > 1)
        spmd(spmdsize)
            if suppress_errors_on_workers
                oldwarnstate0 = warning('off','all');
            end
            % Suppress Matlab error for nearly singular matrix
            oldwarnstate = warning('off','MATLAB:nearlySingularMatrix');
            myname = labindex;
            master = 1; % master process labindex
            mytemp = temps(myname);
            
            if ~isempty(resume)
                S = Sres{myname};
                w = wres{myname};
                M = length(S);
                [prt,pflag] = multpartfunc_w(X,S,w); % Current Tesselation Index
                if pflag
                    error('Tessellation has equivalent centers due to the S and w combination. Change start values.')
                end
                % Log-likelihood
                [llike, ~] = llikefunc(gpgeneral,prt,y,M,Z);
                lprior = get_lprior(n,size(X,2),M,Mmax);
            end
            % Set up a Tessellation Structure
            T = struct('S',S,'M',M,'w',w,'temp',mytemp,...
                'llike',llike,'lprior',lprior);

            % Create independent Random Streams with a seed on each lab
            s = RandStream.create('mrg32k3a','Numstreams',mm,...
                'StreamIndices',myname,'Seed',seed);
            RandStream.setGlobalStream(s);
            if myname == master
                disp('Starting MCMC...')
            end
            for ii=1:(niter+burn)
                % Perform MH step for each parallel region
                [r,T,accepted] = MHstep(y,X,T,gpgeneral,Mmax,n,n_min,sprop_w,Z);
                totals(r) = totals(r) + 1;
                accepts(r) = accepts(r) + accepted;
                naccepted = naccepted + accepted;
                if mod(ii,swapfreq) == 0
                    % Propose a switch of chains and send to all workers
                    if myname == master
                        swapind = zeros(1,2);
                        swapind(1) = randsample(mm-1,1);
                        swapind(2) = swapind(1) + 1;
                        swapind = labBroadcast(master,swapind);
                    else
                        swapind = labBroadcast(master); 
                    end
                    % Send proposed swap to master
                    if myname == swapind(1) && myname ~= master
                        labSend(T,master,1);
                        swaptotal = swaptotal + 1;
                    end
                    if myname == master
                        if myname ~= swapind(1)
                            Tstarswap1 = labReceive(swapind(1),1);
                        else
                            Tstarswap1 = T;
                            swaptotal = swaptotal + 1;
                        end
                    end
                    if myname == swapind(2) && myname ~= master
                        labSend(T,master,2);
                        swaptotal = swaptotal + 1;
                    end
                    if myname == master
                        if myname ~= swapind(2)
                            Tstarswap2 = labReceive(swapind(2),2);
                        else
                            Tstarswap2 = T;
                            swaptotal = swaptotal + 1;
                        end
                    end
                    if myname == master
                        swaptotal_global = swaptotal_global + 1;
                        lrswap = (Tstarswap2.temp * Tstarswap1.llike + Tstarswap1.lprior) + ...
                            (Tstarswap1.temp * Tstarswap2.llike + Tstarswap2.lprior) - ...
                            (Tstarswap1.temp * Tstarswap1.llike + Tstarswap1.lprior) - ...
                            (Tstarswap2.temp * Tstarswap2.llike + Tstarswap2.lprior);
                        if ~isfinite(lrswap)
                            error('Non-finite swap likelihood ratio.')
                        end
                        if lrswap > log(rand) % Accept
                            swapaccept = 1;
                        else
                            swapaccept = 0;
                        end
                        swapaccepttotal_global = swapaccepttotal_global + swapaccept;
                        swapaccept = labBroadcast(master,swapaccept);
                    else
                        swapaccept = labBroadcast(master);
                    end
                    if swapaccept
                        if myname == master
                            if myname ~= swapind(1)
                                labSend(Tstarswap2,swapind(1))
                            else
                                T = Tstarswap2;
                                swapaccepttotal = swapaccepttotal + 1;
                            end
                            if myname ~= swapind(2)
                                labSend(Tstarswap1,swapind(2))
                            else
                                T = Tstarswap1;
                                swapaccepttotal = swapaccepttotal + 1;
                            end
                        elseif any(myname == swapind)
                            swapaccepttotal = swapaccepttotal + 1;
                            T = labReceive(master);
                        end
                        if any(myname == swapind) % Update temperature on Tree
                            T.temp = mytemp;
                        end
                    end
                end
                if printyes
                    if mod(ii,nprint) == 0
                        disp(['i = ',num2str(ii),', ID = ',num2str(myname),', llike = ',num2str(T.llike),...
                            ', accept = ',num2str(naccepted/ii),...
                            ', swapaccept = ',num2str(swapaccepttotal),'/',num2str(swaptotal),...
                            ', Size = ',num2str(length(T.S)),...
                            ', temp = ',num2str(T.temp)]);
                        if myname == master
                            disp(' ');
                        end
                    end
                end
                if ii > burn
                    % Now Store Values
                    Mpost(ii-burn) = T.M;
                    Spost{ii-burn} = T.S;
                    W(ii-burn,:) = T.w;
                    LLIKE(ii-burn) = T.llike;
                    LPRIOR(ii-burn) = T.lprior;
                end
            end
            % Turn on suppressed warnings
            warning(oldwarnstate);
            if suppress_errors_on_workers
                warning(oldwarnstate0);
                % warning('on','all')
            end
            all_acc_percs = accepts./totals;
            acceptancepercent = naccepted/(niter+burn);
            % Save Output
            % Save output
            output = struct('Mpost',Mpost,...
              'Spost',{Spost},'acceptancepercent',acceptancepercent,...
              'all_acc_percs',all_acc_percs,...
              'W',W,'swap_accept',swapaccepttotal/swaptotal,...
              'n_proposed',totals,...
              'n_accepted',accepts,...
              'llike',LLIKE,'lprior',LPRIOR);
            if saveall
                savenames = 1:mm;
            else
                savenames = master;
            end
            if any(myname == savenames)
                fname = strcat(filepath,'mcmc_id',num2str(myname),'.mat');
                swap_percent_global = swapaccepttotal_global/swaptotal_global;
                savedata(fname,output,swap_percent_global);
            end
            output = [];
        end
        output = [];
        
    elseif(spmdsize == 1)
        T = struct('S',S,'M',M,'w',w,'temp',1,'llike',llike,...
            'lprior',lprior);
        % Set random stream with seed;
        s = RandStream('dsfmt19937','Seed',seed);
        RandStream.setGlobalStream(s);
        
        for ii=1:(niter + burn)
            [r,T,accepted] = MHstep(y,X,T,gpgeneral,Mmax,n,n_min,sprop_w,Z);
            totals(r) = totals(r) + 1;
            accepts(r) = accepts(r) + accepted;
            naccepted = naccepted + accepted;
            if ii > burn
                % Now Store Values
                Mpost(ii-burn) = T.M;
                Spost{ii-burn} = T.S;
                W(ii-burn,:) = T.w;
                LLIKE(ii-burn) = T.llike;
                LPRIOR(ii-burn) = T.lprior;
            end
            if printyes
                percs = accepts./totals;
                if ii <= burn && mod(ii,nprint) == 0
                    disp(['Burn-in ',num2str(ii),'/',num2str(burn),...
                      ', log-lik = ',num2str(T.llike),...
                      ', acceptance = ',num2str(naccepted/ii),...
                      ', birth = ',num2str(percs(1)),...
                      ', death = ',num2str(percs(2)),...
                      ', move = ',num2str(percs(3)),...
                      ', weight = ',num2str(percs(4)),...
                      ', Partitions = ',num2str(T.M)]) % print every 10 iterations
                elseif ii > burn && mod(ii,nprint) == 0
                    disp(['Iteration ',num2str(ii-burn),'/',num2str(niter),...
                      ', log-lik = ',num2str(T.llike),...
                      ', acceptance = ',num2str(naccepted/ii),...
                      ', birth = ',num2str(percs(1)),...
                      ', death = ',num2str(percs(2)),...
                      ', move = ',num2str(percs(3)),...
                      ', weight = ',num2str(percs(4)),...
                      ', Partitions = ',num2str(T.M)]) % print every 10 iterations
                end
            end     
        end
        % Function Return value
        swap_percent = accepts./totals; % It is used in the next line
        output = struct('Mpost',Mpost,...
          'Spost',{Spost},'acceptancepercent',naccepted/(niter + burn),...
          'all_acc_percs',swap_percent,...
          'n_proposed',totals,...
          'n_accepted',accepts,...
          'W',W,...
          'llike',LLIKE,'lprior',LPRIOR);
        % Save data
        fname = strcat(filepath,'mcmc.mat');
        save(fname,'output','swap_percent');
    end   
end % End of Function

function savedata(fname,output,swp_perc)
    save(fname,'output','swp_perc')
end
  
  