function [extreme_tSum,iSE_cluster_exttSum, nC_cand] = cluster_regression_multClust(TC1,TC2,TC3,nE,alphaCut,nC_out)

hVal=nan(1,nE);
pVal=nan(1,nE);
fVal=nan(1,nE);

condNames={'cond1','cond2','cond3'};

iSE_cluster=[]; % indices of start and end of epochs
flagUpdate=1;
countCluster=0;

nS=size(TC1,1);
nC=length(condNames);
for e=1:nE
    
    x_reg=[5, 45, 85]';
    b=nan(nS,2);
    for s=1:nS
        y_reg=[TC1(s,e), TC2(s,e), TC3(s,e)]';
        b(s,:) = regress(y_reg,[ones(nC,1) x_reg]);
    end
    [hVal(e),pVal(e),CI,stats] = ttest(b(:,2),0,'alpha',alphaCut);
    tVal(e)=stats.tstat;
    
    if pVal(e)<alphaCut
        hVal(e)=1;
    else
        hVal(e)=0;
    end 
    
    if hVal(e)==1 
        if e<nE
            if flagUpdate
                iSE_cluster=[iSE_cluster nan(2,1)];
                flagUpdate=0;
                countCluster=countCluster+1;
                iSE_cluster(1,countCluster)=e;
            end
        else
            if countCluster==0
                iSE_cluster(1,1)=e;
                iSE_cluster(2,1)=e;
            else
                if flagUpdate
                    iSE_cluster=[iSE_cluster nan(2,1)];
                    flagUpdate=0;
                    countCluster=countCluster+1;
                    iSE_cluster(1,countCluster)=e;
                    iSE_cluster(2,countCluster)=e;
                else
                    iSE_cluster(2,countCluster)=e;
                end
                
            end
        end
    else
        if flagUpdate==0
            flagUpdate=1;
            if e<nE
                iSE_cluster(2,countCluster)=e-1;
            else
                iSE_cluster(2,countCluster)=e;
            end
        end
    end
end

tSum=zeros(1,nC_out); % if there is no cluster - tSum should be zero
nC=size(iSE_cluster,2); % number of overall clusters
for c=1:nC
    tSum(c)=sum(tVal(iSE_cluster(1,c):iSE_cluster(2,c)));
end

nC_cand=min(nC,nC_out);
extreme_tSum=nan(1,nC_out);
iSE_cluster_exttSum=nan(2,nC_out);

for cou=1:nC_cand
    if nC>0
        [max_fSum,icmax_fSum]=max(tSum);
        [min_fSum,icmin_fSum]=min(tSum);
        if abs(max_fSum)>abs(min_fSum)
            extreme_tSum(cou)=max_fSum;
            icext_fSum=icmax_fSum;
        else
            extreme_tSum(cou)=min_fSum;
            icext_fSum=icmin_fSum;
        end
        iSE_cluster_exttSum(:,cou)=iSE_cluster(:,icext_fSum);
    else
        extreme_tSum(cou)=0;
        iSE_cluster_exttSum(:,cou)=nan(2,1);
    end
    tSum(icext_fSum)=[];
    iSE_cluster(:,icext_fSum)=[];
end

end

