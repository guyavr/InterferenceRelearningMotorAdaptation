function plotMeanTimeCourse_Interference_SepWedges_clusterAna(sha,prot,nC_base_forPlot,clusterAna,col)

names_prot = fieldnames(prot);
for i=1:length(names_prot)
    eval([names_prot{i} '=prot.' names_prot{i} ';']);
end

nCond=length(sha); % nConditions per figure
nC_base=nC_base_forPlot;
nC_clamp=length(iClampC1);
nC_nofbpost=length(iNoFbPostC2);

nC_BC=nC_base+nC_clamp;
nC_tot=nC_BC+nC_nofbpost;

baseC=1:nC_base;
clampC=(nC_base+1):nC_BC;
noFbC=(nC_BC+1):nC_tot;

figSize=[50 100 456 400]; % for 55 cycles
xLim=[0 nC_tot];
yLim=[-5 25];
xTick=[0 nC_base_forPlot 25:20:45 55];
yTick=0:10:60;
FaceAlpha=.3;
c_noFb=[1 1 1]*0.9;
lw=1;
ms=6;
afs=24;
tfs=20;

figure('position',figSize)
hold on
patch('Faces',1:4,'Vertices',[nC_BC+.5 yLim(1); nC_tot yLim(1); nC_tot yLim(2); nC_BC+.5 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')
plot(xLim,[0 0],':k','linewidth',2)
plot([nC_base nC_base]+.5,yLim,'-k','linewidth',2)

for c=1:nCond
    ha=sha{c};
    co=col(c,:);
        
    nS=size(ha,1);
    
    m=nanmean(ha); % mean
    sem_size=nanstd(ha)/sqrt(nS); % standard error
    
    sem=[m-sem_size;m+sem_size];
    
    fill([baseC flip(baseC)]',[sem(1,baseC) flip(sem(2,baseC))]',co,'linestyle','none','facealpha',FaceAlpha);
    fill([clampC flip(clampC)]',[sem(1,clampC) flip(sem(2,clampC))]',co,'linestyle','none','facealpha',FaceAlpha);
    fill([noFbC flip(noFbC)]',[sem(1,noFbC) flip(sem(2,noFbC))]',co,'linestyle','none','facealpha',FaceAlpha);
        
    plot(baseC,m(baseC),'o','markersize',ms,'markerfacecolor',co,'markeredgecolor',co,'linewidth',lw)
    plot(clampC,m(clampC),'o','markersize',ms,'markerfacecolor',co,'markeredgecolor',co,'linewidth',lw)
    plot(noFbC,m(noFbC),'o','markersize',ms,'markerfacecolor',co,'markeredgecolor',co,'linewidth',lw)
    
end

for ses=1:2 % 1- adaptation; 2- aftereffect
    
    if ses==1
        nC_out=clusterAna.adaptation.nC_out;
        iSE_cluster_exttSum=clusterAna.adaptation.iSE_cluster_exttSum;
        per_outRange=clusterAna.adaptation.per_outRange;
        
        nC_uptoSess=nC_base;
    else
        nC_out=clusterAna.aftereffect.nC_out;
        iSE_cluster_exttSum=clusterAna.aftereffect.iSE_cluster_exttSum;
        per_outRange=clusterAna.aftereffect.per_outRange;
        
        nC_uptoSess=nC_BC;
    end
    
    for clus=1:nC_out
        
        if per_outRange(clus)<0.017 % for 3 comparisons overall (0.05/3)
            plot(nC_uptoSess+[iSE_cluster_exttSum(1,clus)-.4 iSE_cluster_exttSum(2,clus)+.4],-2*[1 1],'-k','linewidth',8)
        end
    end
    
end

% xlabel('Cycle Number (2 Movements)','fontsize',afs)
xlabel(' ','fontsize',afs)
ylabel('Reach angle (deg)','fontsize',afs)
set(gca,'xtick',xTick,'ytick',yTick,'ticklength',[0 0],'xlim',xLim,'ylim',yLim,'fontsize',tfs)

end

