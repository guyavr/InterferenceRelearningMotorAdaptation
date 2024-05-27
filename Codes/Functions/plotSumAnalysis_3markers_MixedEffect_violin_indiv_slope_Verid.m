function plotSumAnalysis_3markers_MixedEffect_violin_indiv_slope_Verid(a,col,yLabel,yLim,yTick)

as=a.indiv;
nS=size(as,1);
nC=size(as,2);

m=a.m;
se=a.seIntrvl;

inter = a.MixedEffect.fixed.inter;
slope = a.MixedEffect.fixed.slope;

range_as=yLim(1,2)-yLim(1,1);
yTick_main=yTick{1};
figSize=[50 100 300 400];
xLim=[-10, 100];

xTick=[5 45 85];
bw=0.55;
lw=1;
tfs=20;
txtfs=14;

ldist=0.12*range_as; % line distances
txtdist=0.054*range_as; % line distances
nP=3; % number of pVal
l_loc=yLim(1,2)-(1:nP)*ldist;
ytxt_loc=l_loc+txtdist;

% plot means and var of each cond
figure('position',figSize)
hold on

for c=1:nC
    co=col(c,:);  
    plot([1 1]*xTick(c),[se(1,c) se(2,c)],'-','color',[1 1 1]*0.2,'linewidth',3);
    plot(xTick(c),m(c),'o','linewidth',lw, 'markersize', 12, 'markeredgecolor',co,'markerfacecolor',co);
end

x1=xLim(1)+10; x2=xLim(2)-10;
plot([x1,x2],inter+slope*[x1 x2],'-','Color','k','LineWidth',4)
plot(xLim,[0 0],':k','linewidth',2)

set(gca,'xtick',xTick,'ytick',yTick_main,'fontsize',tfs)
ax=gca;
ylabel(yLabel)
xlim(xLim)
ylim(yLim(1,:))

% plotting the comparison between conditions through indiv differences
xshift=0.05;
xloc=0.07.*rand(nS,1) +xshift;

xBar=1.08;
bwDiff=0.09;

slope_indiv = a.MixedEffect.random.slope';
slope_CI = double([a.MixedEffect.fixed.stats(2,7) a.MixedEffect.fixed.stats(2,8)]);

% Individual+violin
xLim_diff_vio= [0.7 1.2];
yLim_diff_vio=yLim(3,:);
yTick_diff_vio=yTick{3};

co=col(2,:);

figure('position',[50 100 180 300])
hold on
plot(xLim_diff_vio,[0 0],':k','linewidth',2)

scatter(0.9-xloc,slope_indiv,50,'o','markeredgecolor','none','markerfacecolor',co,'markerfacealpha',.5);

bandwidth = diff(yLim_diff_vio)*0.0375;

xvalues = linspace( prctile(slope_indiv,1), prctile(slope_indiv,99), 100 );
[f,xi] = ksdensity(slope_indiv,xvalues,'Bandwidth',bandwidth,'BoundaryCorrection','reflection');
patch( 0.9 +[f,zeros(1,numel(xi),1),0]*0.009,[xi,fliplr(xi),xi(1)],co)

% difference bar
bar(xBar,slope,'barwidth',bwDiff,'linewidth',2, 'edgecolor','k','facecolor',co,'ShowBaseLine','off');
plot(xBar*[1 1],[slope_CI(1) slope_CI(2)],'-k','linewidth',2)

set(gca,'xtick',[],'xticklabel',[],'ytick',yTick_diff_vio,'fontsize',tfs,'tickdir','out')
ax=gca;
ax.XAxis.Visible = 'off';
xlim(xLim_diff_vio);
ylim(yLim_diff_vio);
