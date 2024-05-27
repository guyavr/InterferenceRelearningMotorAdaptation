function [extreme_tSum,iSE_cluster_exttSum, nC_cand] = cluster_ttest_multClust(TC1,TC2,nE,alphaCut,nC_out)

hVal=nan(1,nE);
pVal=nan(1,nE);
tVal=nan(1,nE);

iSE_cluster=[]; % indices of start and end of epochs
flagUpdate=1;
countCluster=0;
for e=1:nE
    [hVal(e),pVal(e),~,stats] = ttest(TC2(:,e),TC1(:,e),'alpha',alphaCut);
    tVal(e)=stats.tstat;
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
        [max_tSum,icmax_tSum]=max(tSum);
        [min_tSum,icmin_tSum]=min(tSum);
        if abs(max_tSum)>abs(min_tSum)
            extreme_tSum(cou)=max_tSum;
            icext_tSum=icmax_tSum;
        else
            extreme_tSum(cou)=min_tSum;
            icext_tSum=icmin_tSum;
        end
        iSE_cluster_exttSum(:,cou)=iSE_cluster(:,icext_tSum);
    else
        extreme_tSum(cou)=0;
        iSE_cluster_exttSum(:,cou)=nan(2,1);
    end
    tSum(icext_tSum)=[];
    iSE_cluster(:,icext_tSum)=[];
end

end

