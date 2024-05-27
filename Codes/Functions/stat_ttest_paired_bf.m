function ana = stat_ttest_paired_bf(a)
% The function returns the results of a ttest analysis and the BF10 [base
% factor in favor of H1 (bf10) and in favor of H0 (bf01)].
% a- the data of all participants

a_s=a.indiv;
nS=size(a_s,1);
nC=size(a_s,2);

% bootstrap for non-parametric confidence interval of the mean
nboot=1000;
bootfun=@(x)nanmean(x);

lilh=nan(1,nC);
lilp=nan(1,nC);
jbh=nan(1,nC);
jbp=nan(1,nC);

if nS>=4
    for c=1:nC
        [lilh(c),lilp(c)] = lillietest(a_s(:,c));
        [jbh(c),jbp(c)] = jbtest(a_s(:,c));
    end
end
            
m=nanmean(a_s);
se=nanstd(a_s)/sqrt(nS);
seIntrvl=[m-se;m+se];
ci=1.96*se;
ciRange=[m-ci;m+ci];

% diff scores
a_s_diff=a_s(:,2)-a_s(:,1);
m_diff=nanmean(a_s_diff);
se_diff=nanstd(a_s_diff)/sqrt(nS);
seIntrvl_diff=[m_diff-se_diff;m_diff+se_diff];

[lilh_diff,lilp_diff] = lillietest(a_s_diff);

y=a_s;
[h,p,ci_diff,stats]=ttest(y(:,2),y(:,1));

% based on a Matlab package for Bayes Factor statistical analysis.
% See ..Google Drive\Matlab_FileExchange\bayesFactor-master\bayesFactor-master
[bf10,~] = bf.ttest(y(:,2),y(:,1));

% Cohen's d effect size
cohend = computeCohen_d(y(:,2), y(:,1), 'paired');

ana.indiv=a_s;
ana.m=m;
ana.ci=ciRange;
ana.se=se;
ana.seIntrvl=seIntrvl;
ana.lilh=lilh;
ana.lilp=lilp;
ana.jbh=jbh;
ana.jbp=jbp;

ana.indiv_diff=a_s_diff;
ana.m_diff=m_diff;
ana.ci_diff=ci_diff;
ana.se_diff=se_diff;
ana.seIntrvl_diff=seIntrvl_diff;
ana.lilh_diff=lilh_diff;
ana.lilp_diff=lilp_diff;

ana.ttest.h=h;
ana.ttest.p=p;
ana.ttest.stats=stats;

ana.bf10=bf10;
ana.bf01=1/bf10;

ana.cohend=cohend;

end

