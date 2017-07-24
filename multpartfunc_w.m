% MULTPARTFUNC_W: Assign partition labels to observations (rows) in X
%                 using a weighted Voronoi tessellation.
%
% prt = multpartfunc_w(X,S,w)
%
% X is a nxp design matrix
% S is a Mx1 vector of indices of X which are the tesselation centers
% w is a px1 vector of weights (already normalized such that sum(w.^2) = 1)
%
% prt is a vector of same length as size(X,1) containing the partition 
%   element label for each observation in X.  The element labels are an
%   index on S.

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

function prt = multpartfunc_w(X,S,w)
    n = size(X,1);
    M = size(S,1);
    prt = zeros(n,1);
    for ii=1:n
      dists = zeros(M,1);
      for jj=1:M
        % dists(jj) = norm( X(ii,:) - X(S(jj),:) );
        dists(jj) = ((X(ii,:) - X(S(jj),:)).^2) * (w.^2);
      end
      [~,minind] = min(dists);
      prt(ii) = minind;
    end
    
    if(max(prt) ~= M)
        S
        X(S,:)
        [max(prt), M]
        error('Partitions & M Do not agree! Most likely duplicate values in design matrix.')
    end
end