function plotMeanTimeCourse_allCyc_Interference_AntiCLongV_clusterAna(sha_bothG,prot,col_bothG,clusterAna)

% plot the figure based on the protocol of the test group
names_prot = fieldnames(prot);
for i=1:length(names_prot)
    eval([names_prot{i} '=prot.' names_prot{i} ';']);
end

nG=length(sha_bothG); % nConditions per figure

figSize=[50 100 1173 400];

xLim=[0 nC];
yLim=[-5 25];
xTick=[0 iNoFbC1(end) iFbC1(end) iClampC1(end) iNoFbPostC1(end) iWashC(end) iClampC2(end) iNoFbPostC2(end) ];
yTick=0:10:60;
FaceAlpha=.3;
c_noFb=[1 1 1]*0.9;
lw=1;
tfs=20;

figure('position',figSize)
hold on
patch('Faces',1:4,'Vertices',[0 yLim(1); iNoFbC1(end)+.5 yLim(1); iNoFbC1(end)+.5 yLim(2); 0 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')
patch('Faces',1:4,'Vertices',[iNoFbPostC1(1)-.5 yLim(1); iNoFbPostC1(end)+.5 yLim(1); iNoFbPostC1(end)+.5 yLim(2); iNoFbPostC1(1)-.5 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')
patch('Faces',1:4,'Vertices',[iNoFbPostC2(1)-.5 yLim(1); iNoFbPostC2(end)+.5 yLim(1); iNoFbPostC2(end)+.5 yLim(2); iNoFbPostC2(1)-.5 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')
plot(xLim,[0 0],':k','linewidth',2)

plot([iFbC1(end) iFbC1(end)]+.5,yLim,'-k','linewidth',2)
plot([iWashC(end) iWashC(end)]+.5,yLim,'-k','linewidth',2)

for g=1:nG
    ha=sha_bothG{g};
    co=col_bothG{g};
    
    nS=size(ha,1);
    
    m=nanmean(ha); % mean
    sem_size=nanstd(ha)/sqrt(nS); % standard error
    
    sem=[m-sem_size;m+sem_size];
    
    fill([iNoFbC1 flip(iNoFbC1)]',[sem(1,iNoFbC1) flip(sem(2,iNoFbC1))]',co(1,:),'linestyle','none','facealpha',FaceAlpha);
    fill([iFbC1 flip(iFbC1)]',[sem(1,iFbC1) flip(sem(2,iFbC1))]',co(1,:),'linestyle','none','facealpha',FaceAlpha);
    fill([iClampC1 flip(iClampC1)]',[sem(1,iClampC1) flip(sem(2,iClampC1))]',co(1,:),'linestyle','none','facealpha',FaceAlpha);
    fill([iNoFbPostC1 flip(iNoFbPostC1)]',[sem(1,iNoFbPostC1) flip(sem(2,iNoFbPostC1))]',co(1,:),'linestyle','none','facealpha',FaceAlpha);
    fill([iWashC flip(iWashC)]',[sem(1,iWashC) flip(sem(2,iWashC))]',co(2,:),'linestyle','none','facealpha',FaceAlpha);
    fill([iClampC2 flip(iClampC2)]',[sem(1,iClampC2) flip(sem(2,iClampC2))]',co(2,:),'linestyle','none','facealpha',FaceAlpha);
    fill([iNoFbPostC2 flip(iNoFbPostC2)]',[sem(1,iNoFbPostC2) flip(sem(2,iNoFbPostC2))]',co(2,:),'linestyle','none','facealpha',FaceAlpha);
    
    plot(iNoFbC1,m(iNoFbC1),'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:),'linewidth',lw)
    plot(iFbC1,m(iFbC1),'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:),'linewidth',lw)
    plot(iClampC1,m(iClampC1),'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:),'linewidth',lw)
    plot(iNoFbPostC1,m(iNoFbPostC1),'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:),'linewidth',lw)
    plot(iWashC,m(iWashC),'o','markerfacecolor',co(2,:),'markeredgecolor',co(2,:),'linewidth',lw)
    plot(iClampC2,m(iClampC2),'o','markerfacecolor',co(2,:),'markeredgecolor',co(2,:),'linewidth',lw)
    plot(iNoFbPostC2,m(iNoFbPostC2),'o','markerfacecolor',co(2,:),'markeredgecolor',co(2,:),'linewidth',lw)
    
end

for ses=1:7 % 1- no fb; 2- baseline; 3- adaptation 1; 4- aftereffect 1; 5- washout; 6- adaptation 2; 7- aftereffect 2
    
    xTick=[0 iNoFbC1(end) iFbC1(end) iClampC1(end) iNoFbPostC1(end) iWashC(end) iClampC2(end) iNoFbPostC2(end) ];

    if ses==1
        nC_out=clusterAna.noFb.nC_out;
        iSE_cluster_exttSum=clusterAna.noFb.iSE_cluster_exttSum;
        per_outRange=clusterAna.noFb.per_outRange;
        
        nC_uptoSess=0;
    elseif ses==2
        nC_out=clusterAna.baseline.nC_out;
        iSE_cluster_exttSum=clusterAna.baseline.iSE_cluster_exttSum;
        per_outRange=clusterAna.baseline.per_outRange;
        
        nC_uptoSess=iNoFbC1(end);
    elseif ses==3
        nC_out=clusterAna.adaptation1.nC_out;
        iSE_cluster_exttSum=clusterAna.adaptation1.iSE_cluster_exttSum;
        per_outRange=clusterAna.adaptation1.per_outRange;
        
        nC_uptoSess=iFbC1(end);
    elseif ses==4
        nC_out=clusterAna.aftereffect1.nC_out;
        iSE_cluster_exttSum=clusterAna.aftereffect1.iSE_cluster_exttSum;
        per_outRange=clusterAna.aftereffect1.per_outRange;
        
        nC_uptoSess=iClampC1(end);
    elseif ses==5
        nC_out=clusterAna.washout.nC_out;
        iSE_cluster_exttSum=clusterAna.washout.iSE_cluster_exttSum;
        per_outRange=clusterAna.washout.per_outRange;
        
        nC_uptoSess=iNoFbPostC1(end);
    elseif ses==6
        nC_out=clusterAna.adaptation2.nC_out;
        iSE_cluster_exttSum=clusterAna.adaptation2.iSE_cluster_exttSum;
        per_outRange=clusterAna.adaptation2.per_outRange;
        
        nC_uptoSess=iWashC(end);
    elseif ses==7
        nC_out=clusterAna.aftereffect2.nC_out;
        iSE_cluster_exttSum=clusterAna.aftereffect2.iSE_cluster_exttSum;
        per_outRange=clusterAna.aftereffect2.per_outRange;
        
        nC_uptoSess=iClampC2(end);
    end
    
    for clus=1:nC_out
        
        if per_outRange(clus)<0.017%0.05 % since have 3 conditions overall (0.05/3)
            plot(nC_uptoSess+[iSE_cluster_exttSum(1,clus)-.4 iSE_cluster_exttSum(2,clus)+.4],-2.5*[1 1],'-k','linewidth',8)
        end
    end
    
end

xlabel('Cycle number (2 Movements)')
ylabel('Reach angle (deg)')
set(gca,'xtick',xTick,'ytick',yTick,'ticklength',[0 0],'xlim',xLim,'ylim',yLim,'fontsize',tfs)

end

