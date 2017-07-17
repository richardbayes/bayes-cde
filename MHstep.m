function [r,T,accepted] = MHstep(y,X,T,gpgeneral,Mmax,n,n_min,precision,Z)
    p = size(X,2);
    accepted = 0;
    % Determine which step to do 
    if T.M > 1 && T.M < Mmax
        r = randsample(4,1);
        birthrev = 0;
        deathrev = 0;
    elseif T.M == 1
        r = randsample([1 3 4],1);
    elseif T.M == Mmax
        r = randsample([2 3 4],1);
    else
        error('A case that was not expected occurred in MCMC')
    end
    if r == 1 % Birth
        % Choose an index not already in S
        ind = ismember(1:n,T.S);
        IND = 1:n;
        IND = IND(ind == 0);
        Sprop = [T.S; randsample(IND,1)]; % Proposed Tesselation Center
        Mprop = T.M + 1;
        % proposed partition
        prtprop = multpartfunc_w(X,Sprop,T.w);
        % Make sure we have a minimum number of obs in each partition
        tab = tabulate(prtprop);
        if(all(tab(:,2) > n_min)) % Proceed only if we have enough obs in each partition
            % Calculate log likelihood
            [llikeprop, ~] = llikefunc(gpgeneral,prtprop,y,Mprop,Z);
            lpriorprop = get_lprior(n,p,Mprop,Mmax);
            % Corrections on boundary cases to the log MH ratio
            if T.M > 1 && T.M < (Mmax - 1)
                birthrev = 0;
            elseif T.M == 1
                birthrev = log(3/4);
            elseif T.M == (Mmax - 1)
                birthrev = log(4/3);
            else
                error('Unexpected case in Birth Step')
            end
            % log-MH ratio
            lr = T.temp * (llikeprop - T.llike) + lpriorprop - T.lprior + birthrev;
            if lr > log(rand(1))
                T.llike = llikeprop;
                T.S = Sprop;
                T.M = Mprop;
                T.lprior = lpriorprop;
                % prt = prtprop;
                accepted = 1;
                % GP = GPprop;
            end
        end
    elseif r == 2 % Death Step
        % Randomly delete index from existing centers
        Sprop = T.S;
        Sprop(randsample(size(T.S,1),1)) = [];
        Mprop = T.M - 1;
        % proposed partition
        prtprop = multpartfunc_w(X,Sprop,T.w);

        % Calculate log likelihood
        [llikeprop, ~] = llikefunc(gpgeneral,prtprop,y,Mprop,Z);
        lpriorprop = get_lprior(n,p,Mprop,Mmax);
        % Corrections on boundary cases to the log MH ratio
        if T.M > 2 && T.M < Mmax
            deathrev = 0;
        elseif T.M == 2
            deathrev = log(4/3);
        elseif T.M == Mmax
            deathrev = log(3/4);
        else
            error('Unexpected case in Death Step')
        end
        % log-MH ratio
        lr = T.temp * (llikeprop - T.llike) + lpriorprop - T.lprior + deathrev;
        if lr > log(rand(1))
            T.llike = llikeprop;
            T.lprior = lpriorprop;
            T.S = Sprop;
            T.M = Mprop;
            % prt = prtprop;
            accepted = 1;
            % GP = GPprop;
        end
    elseif r == 3 % Move Step
        mind = randsample(T.M,1);
        Sprop = T.S;
        Sprop(mind) = []; % Delete one index
        % Add in one index from ones not previously in S
        ind = ismember(1:n,T.S);
        IND = 1:n;
        IND = IND(ind == 0);
        Sprop = [Sprop; randsample(IND,1)];

        % Proposed M
        Mprop = T.M;

        % proposed partition
        prtprop = multpartfunc_w(X,Sprop,T.w);
        % Make sure we have a minimum number of obs in each partition
        tab = tabulate(prtprop);
        if(all(tab(:,2) > n_min)) % Proceed only if we have enough obs in each partition
            % Calculate log likelihood
            [llikeprop, ~] = llikefunc(gpgeneral,prtprop,y,Mprop,Z);
            lpriorprop = get_lprior(n,p,Mprop,Mmax);
            % No boundary cases to correct for with a move step

            % log-MH ratio
            lr = T.temp * (llikeprop - T.llike) + lpriorprop - T.lprior;
            if lr > log(rand(1))
                T.llike = llikeprop;
                T.S = Sprop;
                T.M = Mprop;
                T.lprior = lpriorprop;
                % prt = prtprop;
                accepted = 1;
                % GP = GPprop;
            end
        end
    elseif r == 4 % Change weights step
        % wtotal = wtotal + 1;

        % Propose a new weight using a Dirichlet distribution as the
        %   proposal distribution
        % Precision parameter (higher precision increases acceptance rates)
        alpha = precision*T.w.^2;

        wprop = sqrt(drchrnd(alpha',1)');
        alphaprop = precision*wprop.^2;

        prtprop = multpartfunc_w(X,T.S,wprop);
        % Make sure we have a minimum number of obs in each partition
        tab = tabulate(prtprop);
        if(all(tab(:,2) > n_min)) % Proceed only if we have enough obs in each partition
            % Calculate log likelihood
            [llikeprop, ~] = llikefunc(gpgeneral,prtprop,y,T.M,Z);
            lpriorprop = get_lprior(n,p,T.M,Mmax);
            % log-MH ratio
            % Uniform prior means we don't have anything but the likelihood 
            %   to determine acceptance since we have uniform priors on
            %   everything
            % Add in the dirichlet penalty since it is not "symmetric"
            pen = drchpdf(T.w'.^2,alphaprop',1) - drchpdf(wprop'.^2,alpha',1);
            lr = T.temp * (llikeprop - T.llike) + lpriorprop - T.lprior + pen;
            % waccept = 0;
            if lr > log(rand(1))
                T.llike = llikeprop;
                T.w = wprop;
                T.lprior = lpriorprop;
                %prt = prtprop;
                accepted = 1;
                % waccept = 1;
                %GP = GPprop;
            end
        end
    end
end

