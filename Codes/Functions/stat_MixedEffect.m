function ana = stat_MixedEffect(a,cond)
% a- the data of all participants

a_s=a;
nS=size(a_s,1);
nC=size(a_s,2);

% bootstrap for non-parametric confidence interval of the mean
nboot=1000;
bootfun=@(x)nanmean(x);

lilh=nan(1,nC);
lilp=nan(1,nC);

if nS>=4
    for c=1:nC
        [lilh(c),lilp(c)] = lillietest(a_s(:,c));
    end
end
            
m=nanmean(a_s);
ci=bootci(nboot,bootfun,a_s);
se=nanstd(a_s)/sqrt(nS);
seIntrvl=[m-se;m+se];

% LME
y=[a_s(:,1);a_s(:,2);a_s(:,3)];
nV=[cond(1)*ones(nS,1);cond(2)*ones(nS,1);cond(3)*ones(nS,1)];
subj=repmat([1:nS]',length(cond),1);
tbl = table(y,nV,subj,'VariableNames',{'HA','nVerid','Subj'});
lme = fitlme(tbl,'HA ~ nVerid + (nVerid|Subj)'); % estimate both fixed and ramdom effect
[beta,betanames,stats_fixed] = fixedEffects(lme);
[B,Bnames,stats_random]= randomEffects(lme);
B_subj=reshape(B,2,[]); % first raw- intercept, second- slope

writetable(tbl,'Attenuation_Verid_aftereffect_MixedEffect.csv','Delimiter',',')
inter_fixed=beta(1);
slope_fixed=beta(2);

inter_random=beta(1)+ B_subj(1,:);
slope_random=beta(2)+ B_subj(2,:);

ana.indiv=a_s;
ana.m=m;
ana.ci=ci;
ana.se=se;
ana.seIntrvl=seIntrvl;
ana.lilh=lilh;
ana.lilp=lilp;

ana.MixedEffect.fixed.inter=inter_fixed;
ana.MixedEffect.fixed.slope=slope_fixed;
ana.MixedEffect.fixed.stats=stats_fixed;

ana.MixedEffect.random.inter=inter_random;
ana.MixedEffect.random.slope=slope_random;
ana.MixedEffect.random.stats=stats_random;

end

