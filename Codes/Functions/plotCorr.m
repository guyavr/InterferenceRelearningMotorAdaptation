function plotCorr(a,col,axes_label,axes_lim,axes_tick)

a_s=a.indiv;
xLim_corr=axes_lim(1,:);
yLim_corr=axes_lim(2,:);
[rho,pval] = corr(a_s);
corr_r = rho(1,2);
corr_p = pval(1,2);

% Plot correlation
figure('position',[50 100 400 400])
hold on
plot(xLim_corr,[0 0],':k','linewidth',2)
plot([0 0],yLim_corr,':k','linewidth',2)
pc=scatter(a_s(:,1),a_s(:,2),'o','markeredgecolor','none','markerfacecolor',col,'sizedata',120);
alpha 0.7
set(gca,'xlim',xLim_corr,'ylim',yLim_corr,'xtick',axes_tick{1},'ytick',axes_tick{2},'fontsize',20)
xlabel(axes_label{1},'fontsize',20)
ylabel(axes_label{2},'fontsize',20)
axis square

end

