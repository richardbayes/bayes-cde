load('sim6_BPM.mat')
% Look at partition with maximum posterior probability
[~,I] = max(chainstats.LL_store);
% MAP Voronoi centers
centers{I}

% Plot likelihood
plot(chainstats.LL_store)
% Plot number of tessellations
plot(chainstats.k_store)

% Get posteriors with credible intervals
ngrid = 100;
xgrid = linspace(-3,10,ngrid);
nmcmc = size(preds,1);
postdens = [];
for jj=1:size(datpred,1) % For every prediction spot
    tmp = zeros(nmcmc,ngrid);
    for ii=1:nmcmc
      tmp(ii,:) = normpdf(xgrid,themeans(ii,jj),sqrt(sigma2(ii)));
    end
    qtls = quantile(tmp,[.025,.975])';
    pmean = mean(tmp);
    postdens(jj).lb = qtls(:,1);
    postdens(jj).ub = qtls(:,2);
    postdens(jj).m = pmean';
end
% Now plot 
figure()
figlocs = [3 4 1 2];
for jj=1:size(datpred,1)
    subplot(2,2,figlocs(jj));
    plot(xgrid,postdens(jj).m,'--k','LineWidth',2)
    hold on
        plot(xgrid,postdens(jj).lb,'--k')
        plot(xgrid,postdens(jj).ub,'--k')
    hold off
    title(strcat('X_1=',num2str(datpred(jj,2)),...
                 ', X_2=',num2str(datpred(jj,3))));
    xlabel('Y');
    ylabel('density');
    % Plot true function...
    hold on
        if datpred(jj,2) > datpred(jj,3) && datpred(jj,3) < .75
            plot(xgrid,gampdf(xgrid,10,.5),'r')
        elseif datpred(jj,2) < datpred(jj,3) && datpred(jj,2) < .75
            bifun = @(x) .5*normpdf(x,1,1) + .5*normpdf(x,5,1);
            plot(xgrid,bifun(xgrid),'r')
        else
            plot(xgrid,normpdf(xgrid,1,sqrt(.5)),'r');
        end
    hold off
end
%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 6 6];
%saveas(fig,'../../figs/sim6BPM.jpg')
%print('../../figs/sim6BPM','-depsc')
% Save data so the graphic can be made in R (Matlab graphic bugs)
rdat = [xgrid', postdens(1).m, postdens(2).m, postdens(3).m, postdens(4).m, ...
    postdens(1).lb, postdens(1).ub,...
    postdens(2).lb, postdens(2).ub,...
    postdens(3).lb, postdens(3).ub,...
    postdens(4).lb, postdens(4).ub];
csvwrite('Rfigs/data/section4_2_BPM.csv', rdat);