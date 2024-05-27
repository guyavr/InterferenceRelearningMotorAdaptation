function plotExpDesign_Interference_AntiCLongV_LongBaseline(prot,col)
% prot - protocol
% col - colors

names_prot = fieldnames(prot);
for i=1:length(names_prot)
    eval([names_prot{i} '=prot.' names_prot{i} ';']);
end

xLim=[0 nC];
yLim=[-15 15];
c_noFb=[1 1 1]*0.9;
tfs=20;

figSize=[50 100 1173 110];

figure('position',figSize)
hold on

co=mean(col);

patch('Faces',1:4,'Vertices',[0 yLim(1); iNoFbC1(end)+.5 yLim(1); iNoFbC1(end)+.5 yLim(2); 0 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')
patch('Faces',1:4,'Vertices',[iNoFbPostC1(1)-.5 yLim(1); iNoFbPostC1(end)+.5 yLim(1); iNoFbPostC1(end)+.5 yLim(2); iNoFbPostC1(1)-.5 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')
patch('Faces',1:4,'Vertices',[iNoFbPostC2(1)-.5 yLim(1); nC yLim(1); nC yLim(2); iNoFbPostC2(1)-.5 yLim(2)],'FaceColor',c_noFb,'EdgeColor','none')

plot([iNoFbC1(end)+.5 iNoFbPostC1(1)-.5],[0 0],':k','linewidth',2)
plot([iNoFbPostC1(end)+.5 iNoFbPostC2(1)-.5],[0 0],':k','linewidth',2)

plot([iFbC1(1)-.5 iFbC1(end)+.5],[0 0],'-','color',co,'linewidth',9)
plot([iClampC1(1)-.5 iClampC1(end)+.5],[0 0],'-','color',co,'linewidth',9)
plot([iWashC(1)-.5 iWashC(end)-.5],[0 0],'-','color',co,'linewidth',9)

plot([iClampC2(1)-.5 iClampC2(end)+.5],[yLim(1) yLim(1)],'-','color',co,'linewidth',9)

set(gca,'xtick',[],'ytick',yLim(1):15:yLim(2),'xlim',xLim,'ylim',yLim,'fontsize',tfs)
h = gca; h.XAxis.Visible = 'off';

end


