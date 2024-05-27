function ana = stat_OneWayANOVA_Repeated_condDiff_AntiCVerid(a,condNames)
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

% ANOVA
y=a_s;
tbl=table(y(:,1),y(:,2),y(:,3),'VariableNames',condNames);
conds = table(categorical(condNames'),'VariableNames',{'Cond'});
rm = fitrm(tbl,[condNames{1} '-' condNames{end} '~1'],'WithinDesign',conds);
[ranovatbl,~,C,~ ]= ranova(rm,'WithinModel','Cond');

pCond=table2array(ranovatbl('(Intercept):Cond','pValue'));

mltcmp=multcompare(rm,'Cond','ComparisonType','bonferroni');

% effect size: partial eta squared
SS_Cond=table2array(ranovatbl('(Intercept):Cond','SumSq'));
SS_CondErr=table2array(ranovatbl('Error(Cond)','SumSq'));
pEtaSq_Cond=SS_Cond/(SS_Cond+SS_CondErr);

% Bayes Factor- not the right way- See JASP analyses
F_Cond=table2array(ranovatbl('(Intercept):Cond','F'));
df1=table2array(ranovatbl('(Intercept):Cond','DF'));
df2=table2array(ranovatbl('Error(Cond)','DF'));

ana.indiv=a_s;
ana.m=m;
ana.ci=ci;
ana.se=se;
ana.seIntrvl=seIntrvl;
ana.lilh=lilh;
ana.lilp=lilp;

ana.OneWay.ranovatbl=ranovatbl;
ana.OneWay.pCond=pCond;
ana.OneWay.mltcmp=mltcmp;

ana.OneWay.pEtaSq_Cond=pEtaSq_Cond;

% diff scores
nComp=3;
for comp=1:nComp
    if comp==1
        a_s_diff=a_s(:,2)-a_s(:,1);
    elseif comp==2
        a_s_diff=a_s(:,3)-a_s(:,1);
    else
        a_s_diff=a_s(:,3)-a_s(:,2);
    end

    m_diff=nanmean(a_s_diff);
    se_diff=nanstd(a_s_diff)/sqrt(nS);
    seIntrvl_diff=[m_diff-se_diff;m_diff+se_diff];
    
    if comp==1
        ana.condDiff.c2_c1.indiv=a_s_diff;
        ana.condDiff.c2_c1.m=m_diff;
        ana.condDiff.c2_c1.ci=table2array([mltcmp(6,6) mltcmp(6,7)]);
        ana.condDiff.c2_c1.se=se_diff;
        ana.condDiff.c2_c1.seIntrvl=seIntrvl_diff;
    elseif comp==2
        ana.condDiff.c3_c1.indiv=a_s_diff;
        ana.condDiff.c3_c1.m=m_diff;
        ana.condDiff.c3_c1.ci=table2array([mltcmp(1,6) mltcmp(1,7)]);
        ana.condDiff.c3_c1.se=se_diff;
        ana.condDiff.c3_c1.seIntrvl=seIntrvl_diff;
    else
        ana.condDiff.c3_c2.indiv=a_s_diff;
        ana.condDiff.c3_c2.m=m_diff;
        ana.condDiff.c3_c2.ci=table2array([mltcmp(2,6) mltcmp(2,7)]);
        ana.condDiff.c3_c2.se=se_diff;
        ana.condDiff.c3_c2.seIntrvl=seIntrvl_diff;
    end

end

end

