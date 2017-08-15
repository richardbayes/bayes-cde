function [r,T,accepted] = MHstep2(y,X,T,gpgeneral,Mmax,n,n_min,Z,uprop,sprop_w)
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
        % Choose an existing center to split
        ind = randsample(1:size(T.S,1),1);
        u = normrnd(0,uprop,1,p);
        newS = [T.S(ind,:) - u; T.S(ind,:) + u];
        if all(all(newS > 0) & all(newS < 1)) % Make sure proposal is in the convex hull
            Sprop = T.S;
            Sprop(ind,:) = [];
            Sprop = [Sprop; newS];
            Mprop = T.M + 1;
            % proposed partition
            [prtprop,pflag] = multpartfunc_w2(X,Sprop,T.w);
            if ~pflag
                % Make sure we have a minimum number of obs in each partition
                tab = tabulate(prtprop);
                if(all(tab(:,2) > n_min)) % Proceed only if we have enough obs in each partition
                    % Calculate log likelihood
                    [llikeprop, ~] = llikefunc(gpgeneral,prtprop,y,Mprop,Z);
                    lpriorprop = get_lprior2(n,p,Mprop,Mmax);
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
                    lr = T.temp * (llikeprop - T.llike) + lpriorprop - ...
                        T.lprior + birthrev + log(2/(T.M + 1)) + log(2);
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
            end
        end
    elseif r == 2 % Death Step
        % Randomly merge index from existing centers
        Sprop = T.S;
        ind = randsample(size(T.S,1),2);
        mergedS = (T.S(ind(1),:) + T.S(ind(2),:))./2;
        Sprop(ind,:) = [];
        Sprop = [Sprop; mergedS];
        Mprop = T.M - 1;
        if all(mergedS > 0) && all(mergedS < 1)
            % proposed partition
            [prtprop,pflag] = multpartfunc_w2(X,Sprop,T.w);
            if ~pflag
                tab = tabulate(prtprop);
                if(all(tab(:,2) > n_min))
                    % Calculate log likelihood
                    [llikeprop, ~] = llikefunc(gpgeneral,prtprop,y,Mprop,Z);
                    lpriorprop = get_lprior2(n,p,Mprop,Mmax);
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
                    lr = T.temp * (llikeprop - T.llike) + lpriorprop - T.lprior + ... 
                        deathrev + log((T.M*(T.M-1))/(2*(T.M+1))) - log(2);
                    if lr > log(rand(1))
                        T.llike = llikeprop;
                        T.lprior = lpriorprop;
                        T.S = Sprop;
                        T.M = Mprop;
                        % prt = prtprop;
                        accepted = 1;
                        % GP = GPprop;
                    end
                end
            end
        end
    elseif r == 3 % Move Step
        mind = randsample(T.M,1);
        Sprop = T.S;
        Sprop(mind,:) = []; % Delete one index
        % Randomly sample a place in the convex hull
        Sprop = [Sprop; rand(1,p)];
        % Proposed M
        Mprop = T.M;
        % proposed partition
        [prtprop,pflag] = multpartfunc_w2(X,Sprop,T.w);
        if ~pflag
            % Make sure we have a minimum number of obs in each partition
            tab = tabulate(prtprop);
            if(all(tab(:,2) > n_min)) % Proceed only if we have enough obs in each partition
                % Calculate log likelihood
                [llikeprop, ~] = llikefunc(gpgeneral,prtprop,y,Mprop,Z);
                lpriorprop = get_lprior2(n,p,Mprop,Mmax);
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
        end
    elseif r == 4 % Change weights step
        % Randomly select a weight to change
        ind = randsample(p,1);
        wprop = T.w;
        u = normrnd(0,sprop_w);
        v = ((T.w(ind) + u)/(1 + u) - T.w(ind)) / (T.w(ind) - 1); % value needed to come back
        wprop(ind) = wprop(ind) + u;
        wprop = wprop ./ sum(wprop); % Normalize
        if all(wprop > 0)
            [prtprop,pflag] = multpartfunc_w2(X,T.S,wprop);
            if ~pflag
                % Make sure we have a minimum number of obs in each partition
                tab = tabulate(prtprop);
                if(all(tab(:,2) > n_min)) % Proceed only if we have enough obs in each partition
                    % Calculate log likelihood
                    [llikeprop, ~] = llikefunc(gpgeneral,prtprop,y,T.M,Z);
                    lpriorprop = get_lprior2(n,p,T.M,Mmax);
                    % log-MH ratio
                    % Uniform prior means we don't have anything but the likelihood 
                    %   to determine acceptance since we have uniform priors on
                    %   everything
                    lr = T.temp * (llikeprop - T.llike) + lpriorprop - T.lprior + ...
                            log(normpdf(v,0,sprop_w)) - log(normpdf(u,0,sprop_w));
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
    end
end

