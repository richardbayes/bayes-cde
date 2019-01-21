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

function [r,T,accepted] = MHstep(y,X,T,gpgeneral,Mmax,n,n_min,sprop_w,Z)
% Performs a Metropolis-Hastings step to search the posterior tesselation
%   space in the MCMC algorithm.
%
% Inputs
% y: the dependent variable vector of length n
% X: the covariate design matrix
% T: the current tessellation from which a change is proposed
% gpgeneral: the generic LGP GP object
% Mmax: the maximum number of partition elements allowed
% n: the number of observed data points
% n_min: the minimum number of observations allowed in each partition
%   element.
% sprop_w: the proposal variance for the weights
% Z: the grid that discretizes y
%
% Outputs
% r: a vector of length 4 with three zeros and a one.  The placement of the
%   1 indicates which move was proposed: birth, death, move, change.
% T: the new tessellation if the move is accepted, the old tessellation
%   otherwise.
% accepted: integer. 1 if the proposed move was accepted.

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
        % T.w
        Mprop = T.M + 1;
        % proposed partition
        [prtprop,pflag] = multpartfunc_w(X,Sprop,T.w);
        if ~pflag
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
        end
    elseif r == 2 % Death Step
        % Randomly delete index from existing centers
        Sprop = T.S;
        Sprop(randsample(size(T.S,1),1)) = [];
        %T.w
        Mprop = T.M - 1;
        % proposed partition
        [prtprop,pflag] = multpartfunc_w(X,Sprop,T.w);
        if ~pflag
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
        %T.w
        
        % Proposed M
        Mprop = T.M;

        % proposed partition
        [prtprop,pflag] = multpartfunc_w(X,Sprop,T.w);
        if ~pflag
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
        end
    elseif r == 4 % Change weights step
        % Randomly select a weight to change
        ind = randsample(p,1);
        wprop = T.w;
        u = normrnd(0,sprop_w);
        v = ((T.w(ind) + u)/(1 + u) - T.w(ind)) / (T.w(ind) - 1); % value needed to come back (reversibility)
        wprop(ind) = wprop(ind) + u;
        %T.S
        wprop = wprop ./ sum(wprop); % Normalize
        if all(wprop > 0)
            [prtprop,pflag] = multpartfunc_w(X,T.S,wprop);
            if ~pflag
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

