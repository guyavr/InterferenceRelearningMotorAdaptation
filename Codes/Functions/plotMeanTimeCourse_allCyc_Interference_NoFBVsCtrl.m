function plotMeanTimeCourse_allCyc_Interference_NoFBVsCtrl(sha_bothG,prot_bothG,col_bothG)

% plot the figure based on the protocol of the test group
prot=prot_bothG{1,1};
names_prot = fieldnames(prot);
for i=1:length(names_prot)
    eval([names_prot{i} '=prot.' names_prot{i} ';']);
end

nG=length(sha_bothG); % nConditions per figure
figSize=[50 100 1600 400];
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
patch('Faces',1:4,'Vertices',[iNoFbPostC1(1)-.5 yLim(1); iWashC(end)+.5 yLim(1); iWashC(end)+.5 yLim(2); iNoFbPostC1(1)-.5 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')
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

xlabel('Cycle number (4 Movements)')
ylabel('Reach angle (deg)')
set(gca,'xtick',xTick,'ytick',yTick,'ticklength',[0 0],'xlim',xLim,'ylim',yLim,'fontsize',tfs)

end

