function plotMeanTimeCourse_Washout_GradVsCtrl_clusterAna(sha_bothG,prot_bothG,col_bothG,clusterAna)

% plot the figure based on the protocol of the test group
prot=prot_bothG{1,1};
names_prot = fieldnames(prot);
for i=1:length(names_prot)
    eval([names_prot{i} '=prot.' names_prot{i} ';']);
end

nG=length(sha_bothG); % nConditions per figure
figSize=[50 100 587 400];
xLim=[0 length(iWashC)];
yLim=[-5 25];
xTick=[0 iNoFbC1(end) iFbC1(end) iClampC1(end) iNoFbPostC1(end) iWashC(end) iClampC2(end) iNoFbPostC2(end) ];
yTick=0:10:60;
FaceAlpha=.3;
c_noFb=[1 1 1]*0.9;
lw=1;
tfs=20;

figure('position',figSize)
hold on
plot(xLim,[0 0],':k','linewidth',2)

for g=1:nG
    ha=sha_bothG{g};
    co=col_bothG{g};
    
    nS=size(ha,1);
    
    m=nanmean(ha); % mean
    sem_size=nanstd(ha)/sqrt(nS); % standard error
    
    sem=[m-sem_size;m+sem_size];
    
    fill([iWashC flip(iWashC)]'- iNoFbPostC1(end),[sem(1,iWashC) flip(sem(2,iWashC))]',co(2,:),'linestyle','none','facealpha',FaceAlpha);
    
    plot(iWashC- iNoFbPostC1(end),m(iWashC),'o','markerfacecolor',co(2,:),'markeredgecolor',co(2,:),'linewidth',lw)
    
end

    
xTick=[iNoFbPostC1(end) iWashC(end)];

nC_out=clusterAna.washout.nC_out;
iSE_cluster_exttSum=clusterAna.washout.iSE_cluster_exttSum;
per_outRange=clusterAna.washout.per_outRange;

for clus=1:nC_out

    if per_outRange(clus)<(0.05)
        plot([iSE_cluster_exttSum(1,clus)-.4 iSE_cluster_exttSum(2,clus)+.4],-2.5*[1 1],'-k','linewidth',8)
    end

end


xlabel('Cycle number (4 Movements)')
ylabel('Reach angle (deg)')
set(gca,'xtick',xTick-iNoFbPostC1(end),'xticklabel',xTick,'ytick',yTick,'ticklength',[0 0],'xlim',xLim,'ylim',yLim,'fontsize',tfs)

end

