function plotMeanTimeCourse_allCyc_AttenuationVerid_clusterAna(sha_all,prot,col_all, clusterAna)

% plot the figure based on the protocol of the test group
names_prot = fieldnames(prot);
for i=1:length(names_prot)
    eval([names_prot{i} '=prot.' names_prot{i} ';']);
end

nG=length(sha_all); % nConditions per figure
figSize=[50 100 1100 400]; % for 140 cycles
xLim=[0 nC];
yLim=[-5 25];

xTick=[0+1 iNoFbC1(end)+1 iFb_1targ(end)+1 iFb_2targ(end)-.5 iFbC1(end)+1 iClampC1(end)+.5 iNoFbPostC1(end)];
xTickLabel=[0 iNoFbC1(end) iFb_1targ(end) iFb_2targ(end) iFbC1(end) iClampC1(end) iNoFbPostC1(end)];

yTick=0:10:60;
FaceAlpha=.3;
c_noFb=[1 1 1]*0.9;
lw=1;
tfs=20;

nS=size(sha_all{1, 1},1);

figure('position',figSize)
hold on
patch('Faces',1:4,'Vertices',[0 yLim(1); iNoFbC1(end)+.5 yLim(1); iNoFbC1(end)+.5 yLim(2); 0 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')
patch('Faces',1:4,'Vertices',[iNoFbPostC1(1)-.5 yLim(1); iNoFbPostC1(end)+.5 yLim(1); iNoFbPostC1(end)+.5 yLim(2); iNoFbPostC1(1)-.5 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')
plot(xLim,[0 0],':k','linewidth',2)

plot([iFb_1targ(end) iFb_1targ(end)]+.5,yLim,'--k','linewidth',2)
plot([iFb_2targ(end) iFb_2targ(end)]+.5,yLim,'--k','linewidth',2)
plot([iFbC1(end) iFbC1(end)]+.5,yLim,'-k','linewidth',2)

for g=1:nG
    ha=sha_all{g};
    co=col_all(g,:);

    m=nanmean(ha); % mean
    sem_size=nanstd(ha)/sqrt(nS); % standard error

    sem=[m-sem_size;m+sem_size];

    fill([iNoFbC1 flip(iNoFbC1)]',[sem(1,iNoFbC1) flip(sem(2,iNoFbC1))]',co,'linestyle','none','facealpha',FaceAlpha);
    fill([iFbC1 flip(iFbC1)]',[sem(1,iFbC1) flip(sem(2,iFbC1))]',co,'linestyle','none','facealpha',FaceAlpha);
    fill([iClampC1 flip(iClampC1)]',[sem(1,iClampC1) flip(sem(2,iClampC1))]',co,'linestyle','none','facealpha',FaceAlpha);
    fill([iNoFbPostC1 flip(iNoFbPostC1)]',[sem(1,iNoFbPostC1) flip(sem(2,iNoFbPostC1))]',co,'linestyle','none','facealpha',FaceAlpha);
    
    plot(iNoFbC1,m(iNoFbC1),'o','markerfacecolor',co,'markeredgecolor',co,'linewidth',lw)
    plot(iFb_1targ,m(iFb_1targ),'o','markerfacecolor',co,'markeredgecolor',co,'linewidth',lw)
    plot(iFb_2targ,m(iFb_2targ),'o','markerfacecolor',co,'markeredgecolor',co,'linewidth',lw)
    plot(iFb_alltarg,m(iFb_alltarg),'o','markerfacecolor',co,'markeredgecolor',co,'linewidth',lw)
    plot(iClampC1,m(iClampC1),'o','markerfacecolor',co,'markeredgecolor',co,'linewidth',lw)
    plot(iNoFbPostC1,m(iNoFbPostC1),'o','markerfacecolor',co,'markeredgecolor',co,'linewidth',lw)
    
end

for ses=1:3 % 1- baseline; 2- adaptation; 3- aftereffect
    
    if ses==1
        nC_out=clusterAna.baseline.nC_out;
        iSE_cluster_exttSum=clusterAna.baseline.iSE_cluster_exttSum;
        per_outRange=clusterAna.baseline.per_outRange;
        
        nC_uptoSess=85;

    elseif ses==2
        nC_out=clusterAna.adaptation.nC_out;
        iSE_cluster_exttSum=clusterAna.adaptation.iSE_cluster_exttSum;
        per_outRange=clusterAna.adaptation.per_outRange;
        
        nC_uptoSess=90;
    else
        nC_out=clusterAna.aftereffect.nC_out;
        iSE_cluster_exttSum=clusterAna.aftereffect.iSE_cluster_exttSum;
        per_outRange=clusterAna.aftereffect.per_outRange;
        
        nC_uptoSess=130;
    end
    
    for clus=1:nC_out
        
        if per_outRange(clus)<0.05
            plot(nC_uptoSess+[iSE_cluster_exttSum(1,clus)-.4 iSE_cluster_exttSum(2,clus)+.4],-2*[1 1],'-k','linewidth',8)
        end
    end
end

xlabel('Cycle number (2 Movements)')
ylabel('Reach angle (deg)')
set(gca,'xtick',xTick,'XTickLabel',xTickLabel,'ytick',yTick,'ticklength',[0 0],'xlim',xLim,'ylim',yLim,'fontsize',tfs)

end
