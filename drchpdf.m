% DRCHPDF: evaluate the Dirichlet density.
%
% y = drchpdf(x,a,L) gives the value of the density at x with parameter
%                    vector a.  If L is true, the result is logged.

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

function y = drchpdf(x,a,L)
    t1 = gammaln(sum(a))-sum(gammaln(a));
    t2 = sum((a-1).*log(x));
    y = t1 + t2;
    if ~L
      y = exp(t1 + t2);
    end
end
