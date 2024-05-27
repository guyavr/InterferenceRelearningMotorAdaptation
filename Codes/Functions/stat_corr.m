function ana = stat_corr(a)
% a- the data of all participants

a_s=a.indiv;
nS=size(a_s,1);

% bootstrap for non-parametric confidence interval of the mean
nboot=1000;
bootfun=@(x)nanmean(x);
            
m=nanmean(a_s);
ci=bootci(nboot,bootfun,a_s);
se=nanstd(a_s)/sqrt(nS);
seIntrvl=[m-se;m+se];

y=a_s;
[r,p]=corr(y,'rows','complete');
[corrPear_bf10,~,~] = bf.corr(y(:,1),y(:,2));
    
ana.indiv=a_s;
ana.m=m;
ana.ci=ci;
ana.se=se;
ana.seIntrvl=seIntrvl;

ana.corr.r=r(1,2);
ana.corr.p=p(1,2);
ana.corr.bf10=corrPear_bf10;

end

