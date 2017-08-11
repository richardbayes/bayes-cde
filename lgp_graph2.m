% LGP_GRAPH: plot the results from the output of lgp_partition_model.
%
%     lgp_graph(A,y,X) takes the output A from the lgp_partition_model 
%        function and plots the densities in each partition element
%        from the partition with the highest marginal likelihood.
%        The centers of the partition for which the density is plotted 
%        is listed in the title of each graph. 
%
%        The arguments y and X are the dependent and independent variables
%        originally supplied to lgp_partition_model.
%        
%        If X has one column, a plot of y against X is also produced 
%          with colors indicating different elements of the partition.
%        
%        If X has two columns, the second column of X is plotted against
%          the first column of X, with colors indicating the partition.

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

function lgp_graph2(A,y,X)
    % Find the tessellation with the highest marginal likelihood
    [~,I] = max(A.llike);
    
    % Plot the tesselation if p <= 2
    p = size(X,2);
    % Standardize X for tessellation calculation
    % Xz = zscore(X);
    Xz = X;
    for ii=1:size(X,2)
       Xz(:,ii) = (X(:,ii) - min(X(:,ii))) / (max(X(:,ii)) - min(X(:,ii))); 
    end
    
    % Obtain the partition indices
    prt = multpartfunc_w2(Xz,A.Spost{I},A.W(I,:)');
    if p == 1
        scatter(X,y,30,prt,'filled')
        xlabel('x'); ylabel('y');
    elseif p == 2
        scatter(X(:,1),X(:,2),30,prt,'filled')
        xlabel('x_1'); ylabel('x_2');
    end
    
    % Plot the density in each partition.
    nprt = A.Mpost(I);
    for ii=1:nprt
        figure(ii + 1);
        lgpdens(y(prt == ii));
        title(['Center: ',num2str(A.Spost{I}(ii,:))]);
        ylabel(['Density in Partition ',num2str(ii)])
    end
end