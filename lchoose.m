%     lchoose(n,k) returns the log of the combination n choose k.
% 
%     This file is part of bayes-treed-cde.
% 
%     bayes-treed-cde is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     bayes-treed-cde is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with bayes-treed-cde.  If not, see <http://www.gnu.org/licenses/>.
%
%     Copyright 2017, Richard Payne
function a = lchoose(n,k)
    a = 0;
    for ii=1:n
        a = a + log(ii);
    end
    for ii=1:k
        a = a - log(ii);
    end
    for ii=1:(n-k)
        a = a - log(ii);
    end
end