function plotSumAnalysis_3within_bar_violin_indiv_diff_AntiCLongV(a,col,yLabel,yLim,yTick)

as=a.indiv;
nS=size(as,1);
nC=size(as,2);

m=a.m;
se=a.seIntrvl;
p_t2_g2=a.OneWay.mltcmp.pValue(2);
p_t1_t2=a.OneWay.mltcmp.pValue(4);
p_t1_g2=a.OneWay.mltcmp.pValue(3);

range_as=yLim(1,2)-yLim(1,1);
yTick_main=yTick{1};
figSize=[50 100 300 400];
xLim=[0.5, 3.5];

xTick=1:3;
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
    
    bar(xTick(c),m(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor',co);
    plot([1 1]*xTick(c),[se(1,c) se(2,c)],'-k','linewidth',3);
end

plot(xLim,[0 0],':k','linewidth',2)

plot([xTick(2) xTick(3)],[1 1]*l_loc(3),'-k','linewidth',3)
xtxt_loc=mean([xTick(2) xTick(3)]);
if p_t2_g2<0.05
    if p_t2_g2<0.001
        text(xtxt_loc,ytxt_loc(3),'p<0.001','fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    else
        text(xtxt_loc,ytxt_loc(3),sprintf('p=%.3f',p_t2_g2),'fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    end
else
    text(xtxt_loc,ytxt_loc(3),sprintf('p=%.3f',p_t2_g2),'fontsize',txtfs,'horizontalalignment','center')
end

plot([xTick(1) xTick(2)],[1 1]*l_loc(2),'-k','linewidth',3)
xtxt_loc=mean([xTick(1) xTick(2)]);
if p_t1_t2<0.05
    if p_t1_t2<0.001
        text(xtxt_loc,ytxt_loc(2),'p<0.001','fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    else
        text(xtxt_loc,ytxt_loc(2),sprintf('p=%.3f',p_t1_t2),'fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    end
else
    text(xtxt_loc,ytxt_loc(2),sprintf('p=%.3f',p_t1_t2),'fontsize',txtfs,'horizontalalignment','center')
end

plot([xTick(1) xTick(3)],[1 1]*l_loc(1),'-k','linewidth',3)
xtxt_loc=mean([xTick(1) xTick(3)]);
if p_t1_g2<0.05
    if p_t1_g2<0.001
        text(xtxt_loc,ytxt_loc(1),'p<0.001','fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    else
        text(xtxt_loc,ytxt_loc(1),sprintf('p=%.3f',p_t1_g2),'fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    end
else
    text(xtxt_loc,ytxt_loc(1),sprintf('p=%.3f',p_t1_g2),'fontsize',txtfs,'horizontalalignment','center')
end

set(gca,'xtick',[],'xticklabel',[],'ytick',yTick_main,'fontsize',tfs)
ax=gca;
ylabel(yLabel)
xlim(xLim)
ylim(yLim(1,:))

xshift=0.05;
xloc=0.1.*rand(nS,1) +xshift;

% plotting the comparison between conditions through indiv differences
xshift=0.05;
xloc=0.07.*rand(nS,1) +xshift;

xBar=1.08;
bwDiff=0.09;

nComp=3; % number of pairwise comparisons
for c=1:nComp
    if c==1
        compData=a.condDiff.c2_c1;
        co=col(2,:);
    elseif c==2
        compData=a.condDiff.c3_c1;
        co=col(3,:);
    else
        compData=a.condDiff.c3_c2; 
        co=mean([col(2,:);col(3,:)]); % no need to present this comparison
    end

    m_diff=compData.m;
    v_diff=compData.ci;
    as_diff=compData.indiv;

    % Individual+violin
    xLim_diff_vio= [0.7 1.2];
    yLim_diff_vio=yLim(3,:);
    yTick_diff_vio=yTick{3};

    figure('position',[50 100 180 300])
    hold on
    plot(xLim_diff_vio,[0 0],':k','linewidth',2)

    scatter(0.9-xloc,as_diff,50,'o','markeredgecolor','none','markerfacecolor',co,'markerfacealpha',.5);
    
    bandwidth = diff(yLim_diff_vio)*0.0375;

    xvalues = linspace( prctile(as_diff,1), prctile(as_diff,99), 100 );
    [f,xi] = ksdensity(as_diff,xvalues,'Bandwidth',bandwidth,'BoundaryCorrection','reflection');
    patch( 0.9+[f,zeros(1,numel(xi),1),0],[xi,fliplr(xi),xi(1)],co)
    
    % difference bar
    bar(xBar,m_diff,'barwidth',bwDiff,'linewidth',2, 'edgecolor','k','facecolor',co,'ShowBaseLine','off');
    plot(xBar*[1 1],[v_diff(1) v_diff(2)],'-k','linewidth',2)

    set(gca,'xtick',[],'xticklabel',[],'ytick',yTick_diff_vio,'fontsize',tfs,'tickdir','out')
    ax=gca;
    ax.XAxis.Visible = 'off';
    xlim(xLim_diff_vio);
    ylim(yLim_diff_vio);
end