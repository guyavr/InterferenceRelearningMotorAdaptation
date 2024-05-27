%% Interference - Data Analysis - Exp 4- Anti-Clamp Washout vs Long Veridical FB
clc; clear; close all

addpath('../Functions')

preprocdataDir='../../Data/';

load([preprocdataDir 'AttenuationInterference_AntiClamp_Verid85_trials']);
c_t= [104, 190, 144; 78, 143, 108]/255; % cyan (wedge t: two learning blocks separated by washout)
c_g= [255, 165, 75; 255, 165, 75]/255; % orange (wedge g: one learning block preceded by long veridical feedback) 

nS=max(T.SN);
nT=max(T.TN);
nC=max(T.CN);

% Protocol- hard coded
nC_noFb1=5;
nC_FbC1=5;
nC_ClampC1=40;
nC_NoFbPostC1=10;
nC_Wash=40; % including both anti-clamp and FB baseline
nC_ClampC2=40;
nC_NoFbPostC2=10;

prot=struct;
prot.nC=nC;
prot.iNoFbC1=1:nC_noFb1;
prot.iFbC1=prot.iNoFbC1(end)+(1:nC_FbC1);
prot.iClampC1=prot.iFbC1(end)+(1:nC_ClampC1);
prot.iNoFbPostC1=prot.iClampC1(end)+(1:nC_NoFbPostC1);
prot.iWashC=prot.iNoFbPostC1(end)+(1:nC_Wash);
prot.iClampC2=prot.iWashC(end)+(1:nC_ClampC2);
prot.iNoFbPostC2=prot.iClampC2(end)+(1:nC_NoFbPostC2);

tpc=nT/nC; % trial per cycle

clampCCW=reshape(T.CCW,nT,[])';
clampCCW_s=clampCCW(:,1);

haT_all=T.hand_theta;

% create a matrix of mea (nS X nT)
haT_ccw_cw=reshape(haT_all,nT,[])';
% switch hand angle direction for ccw clamp conditions
haT_cw=haT_ccw_cw;
haT_cw((clampCCW_s==1),:)=-haT_ccw_cw((clampCCW_s==1),:);

% Reaction time
rtT_all=T.RT; % msec
rtT=reshape(rtT_all,nT,[])';

% Movement time
mtT_all=T.MT; % msec
mtT=reshape(mtT_all,nT,[])';

% Search time
stT_all=T.ST;
stT=reshape(stT_all,nT,[])';

% inter trial interval is basically the search time
itiT = stT;

% total trial time- RT, MT and iti
trialT = rtT_all + mtT_all + stT_all;
medianTrialT = median(trialT)/1000; % sec
iqrTrialT = [quantile(trialT,0.25) quantile(trialT,0.75)]/1000; % sec

% rotation
rot_all=T.ri;
rot_all_s=reshape(rot_all,nT,[])';
rot_all_s_ccw=rot_all_s;
rot_all_s_ccw((clampCCW_s==0),:)=-rot_all_s((clampCCW_s==0),:);
% transform to cycles
rot_cyc=nan(nS,nC);

% target angle
ta_all=T.ti;
ta_all_s=reshape(ta_all,nT,[])';

%%
trainWedge=reshape(T.TW,nT,[])'; % wedge with two learning blocks
trainWedge_s=trainWedge(:,1);

ha.indivAllCyc_tw=nan(nS,nC); % individuals' time course: mean hand angle of cycles- training wedge
ha.indivAllCyc_gw=nan(nS,nC); % individuals' time course: mean hand angle of cycles- veridical wedge

nT_out_allSubjRuns=nan(nS,tpc);
nT_org_allSubjRuns=nan(nS,tpc);

for s=1:nS
    % separate ha by targets
    ta_s=ta_all_s(s,:);
    i_nan_ta=find(isnan(ta_s));

    uni_ta=unique(ta_s(isfinite(ta_s))); % output in sorted order

    if trainWedge_s(s)==225
        uni_ta=[uni_ta(3:end) uni_ta(1:2)];
    end

    ha_allTa=nan(tpc,nC);

    for t=1:tpc
        i_ta=find(ta_s==uni_ta(t));
        if length(i_ta)<nC
            ta_s(i_nan_ta(1))=uni_ta(t);
            i_nan_ta(1)=[];
        end
        i_ta_new=find(ta_s==uni_ta(t));

        ha_sepTa_org=haT_cw(s,i_ta_new);
        if length(ha_sepTa_org)<nC
            ha_sepTa_org=[ha_sepTa_org nan(1,nC-length(ha_sepTa_org))];
        end

        [ha_sepTa, nT_org_allSubjRuns(s,t), nT_out_allSubjRuns(s,t)]=removeOutlierTrials_Interference(ha_sepTa_org);  % outlier removal

        ha_allTa(t,:)=ha_sepTa;
    end

    rot_s=rot_all_s_ccw(s,:)';
    rot_s_reshape=reshape(rot_s, tpc,[]);
    rot_cyc(s,:)=rot_s_reshape(1,:);

    ha.indivAllCyc_tw(s,:)=nanmean(ha_allTa(1:2,:));
    ha.indivAllCyc_gw(s,:)=nanmean(ha_allTa(3:4,:));
end

percRemoved=100*sum(sum(nT_out_allSubjRuns))/sum(sum(nT_org_allSubjRuns));

%%

% Separate sessions to Clamp 1 and Clamp 2- and extract measures
nC_base=5;

iBaseClamp1=[prot.iFbC1((end-nC_base+1):end) prot.iClampC1 prot.iNoFbPostC1];
iWashClamp2=[prot.iWashC((end-nC_base+1):end) prot.iClampC2 prot.iNoFbPostC2];

EA_3_7=3:7;

nC_asy=10;
asyCyc=length(iBaseClamp1)-length(prot.iNoFbPostC1)-nC_asy+(1:nC_asy);
AE_all_cyc=asyCyc(end)+(1:length(prot.iNoFbPostC1)); % all no feedback trials

% time course base+clamp
% without baseline subtraction
ha.TC_tw_C1 = ha.indivAllCyc_tw(:,iBaseClamp1);
ha.TC_gw_C1 = ha.indivAllCyc_gw(:,iBaseClamp1);
ha.TC_tw_C2 = ha.indivAllCyc_tw(:,iWashClamp2);
ha.TC_gw_C2 = ha.indivAllCyc_gw(:,iWashClamp2);

% Aftereffect
ha.mAftereffect_all.indiv=[nanmean(ha.TC_tw_C1(:,AE_all_cyc),2),...
    nanmean(ha.TC_gw_C1(:,AE_all_cyc),2),...
    nanmean(ha.TC_tw_C2(:,AE_all_cyc),2),...
    nanmean(ha.TC_gw_C2(:,AE_all_cyc),2)];

%% Statistical Analyses- within groups

% % Hand Angle

% One-way repeated measures ANOVA to enable comparisons across all
% conditions - the gw_c1 is not relevant
condNames={'tw_c1','tw_c2','gw_c2'};

ha.mAftereffect_all = stat_OneWayANOVA_Repeated_condDiff_AntiCVerid(ha.mAftereffect_all.indiv(:,[1 3:4]),condNames);

%% Cluster analysis to identify significant cluster (cycles) differences between the two wedges throughout all cycles of the experiment

newClusterAna_TrainVsGen=0;

if newClusterAna_TrainVsGen

    nCperE=1; % number of cycles per epoch
    
    for bComp=1:7 % 1- no fb; 2- baseline; 3- adaptation 1; 4- aftereffect 1; 5- washout; 6- adaptation 2; 7- aftereffect 2

        if bComp==1
            mC_b1_subj=ha.indivAllCyc_tw(:,prot.iNoFbC1);
            mC_b2_subj=ha.indivAllCyc_gw(:,prot.iNoFbC1);
        elseif bComp==2
            mC_b1_subj=ha.indivAllCyc_tw(:,prot.iFbC1);
            mC_b2_subj=ha.indivAllCyc_gw(:,prot.iFbC1);
        elseif bComp==3
            mC_b1_subj=ha.indivAllCyc_tw(:,prot.iClampC1);
            mC_b2_subj=ha.indivAllCyc_gw(:,prot.iClampC1);
        elseif bComp==4
            mC_b1_subj=ha.indivAllCyc_tw(:,prot.iNoFbPostC1);
            mC_b2_subj=ha.indivAllCyc_gw(:,prot.iNoFbPostC1);
        elseif bComp==5
            mC_b1_subj=ha.indivAllCyc_tw(:,prot.iWashC);
            mC_b2_subj=ha.indivAllCyc_gw(:,prot.iWashC);
        elseif bComp==6
            mC_b1_subj=ha.indivAllCyc_tw(:,prot.iClampC2);
            mC_b2_subj=ha.indivAllCyc_gw(:,prot.iClampC2);
        elseif bComp==7
            mC_b1_subj=ha.indivAllCyc_tw(:,prot.iNoFbPostC2);
            mC_b2_subj=ha.indivAllCyc_gw(:,prot.iNoFbPostC2);
        end

        nS=size(mC_b2_subj,1);
        nC_sess=size(mC_b2_subj,2);
        nE=nC_sess/nCperE;  % number of epochs

        HA_C1=nan(nS,nE);
        HA_C2=nan(nS,nE);

        if nCperE>1
            for e=1:nE
                cE=(e-1)*nCperE+(1:nCperE); % trial number in the epoch
                HA_C1(:,e)=nanmean(mC_b1_subj(:,cE),2);
                HA_C2(:,e)=nanmean(mC_b2_subj(:,cE),2);
            end
        else
            HA_C1=mC_b1_subj;
            HA_C2=mC_b2_subj;
        end

        % Find the tVal of the cluster 
        alphaCut=0.05;
        nC_out=12; % maximum number of clusters to look for
        [extreme_tSum,iSE_cluster_exttSum, nC_cand] = cluster_ttest_multClust(HA_C1,HA_C2,nE,alphaCut,nC_out);

        % pertmutations to the order of the conditions to calculate a
        % distribution of tSum
        nPerm=10000;
        permOpts=randi([0 1], nS,nPerm);
        extreme_tSum_p=nan(nPerm,nC_out);
        iSE_cluster_exttSum_p=cell(nPerm,1);
        for p=1:nPerm
            % Switch between the conditions in the chosen subjects for each
            % permutation
            HA_C1_p=HA_C1;
            HA_C2_p=HA_C2;

            sCondSwitch=find(permOpts(:,p));
            HA_C1_p(sCondSwitch,:)=HA_C2(sCondSwitch,:);
            HA_C2_p(sCondSwitch,:)=HA_C1(sCondSwitch,:);

            [extreme_tSum_p(p,:),iSE_cluster_exttSum_p{p}, ~] = cluster_ttest_multClust(HA_C1_p,HA_C2_p,nE,alphaCut,nC_out);

        end

        per_outRange=nan(1,nC_out);
        for c=1:nC_cand
            
            if extreme_tSum(c)>0
                n_outRange=find( extreme_tSum_p(:,1) > extreme_tSum(c));
            else
                n_outRange=find( extreme_tSum_p(:,1) < extreme_tSum(c));
            end

            if n_outRange
                per_outRange(c)=length(n_outRange)/nPerm;
            else
                per_outRange(c)=0;
            end

        end

        if bComp==1
            clusterAnaResults_TrainVsGen.noFb.nC_out=nC_cand;
            clusterAnaResults_TrainVsGen.noFb.iSE_cluster_exttSum=iSE_cluster_exttSum;
            clusterAnaResults_TrainVsGen.noFb.per_outRange=per_outRange;
            clusterAnaResults_TrainVsGen.noFb.extreme_tSum=extreme_tSum;
            clusterAnaResults_TrainVsGen.noFb.extreme_tSum_p=extreme_tSum_p;
        elseif bComp==2
            clusterAnaResults_TrainVsGen.baseline.nC_out=nC_cand;
            clusterAnaResults_TrainVsGen.baseline.iSE_cluster_exttSum=iSE_cluster_exttSum;
            clusterAnaResults_TrainVsGen.baseline.per_outRange=per_outRange;
            clusterAnaResults_TrainVsGen.baseline.extreme_tSum=extreme_tSum;
            clusterAnaResults_TrainVsGen.baseline.extreme_tSum_p=extreme_tSum_p;
        elseif bComp==3
            clusterAnaResults_TrainVsGen.adaptation1.nC_out=nC_cand;
            clusterAnaResults_TrainVsGen.adaptation1.iSE_cluster_exttSum=iSE_cluster_exttSum;
            clusterAnaResults_TrainVsGen.adaptation1.per_outRange=per_outRange;
            clusterAnaResults_TrainVsGen.adaptation1.extreme_tSum=extreme_tSum;
            clusterAnaResults_TrainVsGen.adaptation1.extreme_tSum_p=extreme_tSum_p;
        elseif bComp==4
            clusterAnaResults_TrainVsGen.aftereffect1.nC_out=nC_cand;
            clusterAnaResults_TrainVsGen.aftereffect1.iSE_cluster_exttSum=iSE_cluster_exttSum;
            clusterAnaResults_TrainVsGen.aftereffect1.per_outRange=per_outRange;
            clusterAnaResults_TrainVsGen.aftereffect1.extreme_tSum=extreme_tSum;
            clusterAnaResults_TrainVsGen.aftereffect1.extreme_tSum_p=extreme_tSum_p;
        elseif bComp==5
            clusterAnaResults_TrainVsGen.washout.nC_out=nC_cand;
            clusterAnaResults_TrainVsGen.washout.iSE_cluster_exttSum=iSE_cluster_exttSum;
            clusterAnaResults_TrainVsGen.washout.per_outRange=per_outRange;
            clusterAnaResults_TrainVsGen.washout.extreme_tSum=extreme_tSum;
            clusterAnaResults_TrainVsGen.washout.extreme_tSum_p=extreme_tSum_p;
        elseif bComp==6
            clusterAnaResults_TrainVsGen.adaptation2.nC_out=nC_cand;
            clusterAnaResults_TrainVsGen.adaptation2.iSE_cluster_exttSum=iSE_cluster_exttSum;
            clusterAnaResults_TrainVsGen.adaptation2.per_outRange=per_outRange;
            clusterAnaResults_TrainVsGen.adaptation2.extreme_tSum=extreme_tSum;
            clusterAnaResults_TrainVsGen.adaptation2.extreme_tSum_p=extreme_tSum_p;
        elseif bComp==7
            clusterAnaResults_TrainVsGen.aftereffect2.nC_out=nC_cand;
            clusterAnaResults_TrainVsGen.aftereffect2.iSE_cluster_exttSum=iSE_cluster_exttSum;
            clusterAnaResults_TrainVsGen.aftereffect2.per_outRange=per_outRange;
            clusterAnaResults_TrainVsGen.aftereffect2.extreme_tSum=extreme_tSum;
            clusterAnaResults_TrainVsGen.aftereffect2.extreme_tSum_p=extreme_tSum_p;
        end

    end

    save('clusterAnaResults_AntiClamp_V85_allCyc','clusterAnaResults_TrainVsGen')

else

    load('clusterAnaResults_AntiClamp_V85_allCyc')

end

%% Cluster analysis to identify significant cluster (cycles) differences between learning functions (clamps)

newClusterAna_clamp=0;

if newClusterAna_clamp
    
    nCperE=1; % number of cycles per epoch
    
    for bComp=1:2 % block comparisions: 1- between two learning block in tw; 2- between 1st learning block in tw and the learning in gw;

        if bComp==1
            mC_b1_subj_all=ha.TC_tw_C1;
            mC_b2_subj_all=ha.TC_tw_C2;
        else
            mC_b1_subj_all=ha.TC_tw_C1;
            mC_b2_subj_all=ha.TC_gw_C2;
        end

        for block=1:3 % 1- baseline; 2- adaptation; 3- aftereffect

            if block==1
                mC_b1_subj=mC_b1_subj_all(:,1:nC_base);
                mC_b2_subj=mC_b2_subj_all(:,1:nC_base);
            elseif block==2
                mC_b1_subj=mC_b1_subj_all(:,nC_base+(1:length(prot.iClampC1)));
                mC_b2_subj=mC_b2_subj_all(:,nC_base+(1:length(prot.iClampC1)));
            else
                mC_b1_subj=mC_b1_subj_all(:,nC_base+length(prot.iClampC1)+1:end);
                mC_b2_subj=mC_b2_subj_all(:,nC_base+length(prot.iClampC1)+1:end);
            end
            
            nS=size(mC_b2_subj,1);
            nC_sess=size(mC_b2_subj,2);
            nE=nC_sess/nCperE;  % number of epochs

            HA_C1=nan(nS,nE);
            HA_C2=nan(nS,nE);

            if nCperE>1
                for e=1:nE
                    cE=(e-1)*nCperE+(1:nCperE); % trial number in the epoch
                    HA_C1(:,e)=nanmean(mC_b1_subj(:,cE),2);
                    HA_C2(:,e)=nanmean(mC_b2_subj(:,cE),2);
                end
            else
                HA_C1=mC_b1_subj;
                HA_C2=mC_b2_subj;
            end

            % Find the tVal of the cluster
            alphaCut=0.05;
            nC_out=12; % maximum number of clusters to look for
            [extreme_tSum,iSE_cluster_exttSum, nC_cand] = cluster_ttest_multClust(HA_C1,HA_C2,nE,alphaCut,nC_out);

            % pertmutations to the order of the conditions to calculate a
            % distribution of tSum
            nPerm=10000;
            permOpts=randi([0 1], nS,nPerm);
            extreme_tSum_p=nan(nPerm,nC_out);
            iSE_cluster_exttSum_p=cell(nPerm,1);
            for p=1:nPerm
                % Switch between the conditions in the chosen subjects for each
                % permutation
                HA_C1_p=HA_C1;
                HA_C2_p=HA_C2;

                sCondSwitch=find(permOpts(:,p));
                HA_C1_p(sCondSwitch,:)=HA_C2(sCondSwitch,:);
                HA_C2_p(sCondSwitch,:)=HA_C1(sCondSwitch,:);

                [extreme_tSum_p(p,:),iSE_cluster_exttSum_p{p}, ~] = cluster_ttest_multClust(HA_C1_p,HA_C2_p,nE,alphaCut,nC_out);

            end

            per_outRange=nan(1,nC_out);
            for c=1:nC_cand
                if extreme_tSum(c)>0
                    n_outRange=find( extreme_tSum_p(:,1) > extreme_tSum(c));
                else
                    n_outRange=find( extreme_tSum_p(:,1) < extreme_tSum(c));
                end

                if n_outRange
                    per_outRange(c)=length(n_outRange)/nPerm;
                else
                    per_outRange(c)=0;
                end
            end

            if block==1
                clusterAnaResults_Clamp.baseline.nC_out=nC_cand;
                clusterAnaResults_Clamp.baseline.iSE_cluster_exttSum=iSE_cluster_exttSum;
                clusterAnaResults_Clamp.baseline.per_outRange=per_outRange;
                clusterAnaResults_Clamp.baseline.extreme_tSum=extreme_tSum;
                clusterAnaResults_Clamp.baseline.extreme_tSum_p=extreme_tSum_p;
            elseif block==2
                clusterAnaResults_Clamp.adaptation.nC_out=nC_cand;
                clusterAnaResults_Clamp.adaptation.iSE_cluster_exttSum=iSE_cluster_exttSum;
                clusterAnaResults_Clamp.adaptation.per_outRange=per_outRange;
                clusterAnaResults_Clamp.adaptation.extreme_tSum=extreme_tSum;
                clusterAnaResults_Clamp.adaptation.extreme_tSum_p=extreme_tSum_p;
            else
                clusterAnaResults_Clamp.aftereffect.nC_out=nC_cand;
                clusterAnaResults_Clamp.aftereffect.iSE_cluster_exttSum=iSE_cluster_exttSum;
                clusterAnaResults_Clamp.aftereffect.per_outRange=per_outRange;
                clusterAnaResults_Clamp.aftereffect.extreme_tSum=extreme_tSum;
                clusterAnaResults_Clamp.aftereffect.extreme_tSum_p=extreme_tSum_p;
            end

        end

        if bComp==1
            clusterAnaResults_Clamp_t1_t2=clusterAnaResults_Clamp;
            save('clusterAnaResults_AntiClamp_V85_Clamp_t1_t2','clusterAnaResults_Clamp_t1_t2')
        else
            save('clusterAnaResults_AntiClamp_V85_Clamp_t1_g2','clusterAnaResults_Clamp_t1_g2')
        end

    end

else
    
    load('clusterAnaResults_AntiClamp_V85_Clamp_t1_t2')
    load('clusterAnaResults_AntiClamp_V85_Clamp_t1_g2')

end
%% Summary time course plots

plotExpDesign_Interference_AntiCLongV_StandardWash(prot,c_t)
plotExpDesign_Interference_AntiCLongV_LongBaseline(prot,c_g)

sha_both={ha.indivAllCyc_tw; ha.indivAllCyc_gw};
col_both={c_t; c_g};

plotMeanTimeCourse_allCyc_Interference_AntiCLongV_clusterAna(sha_both,prot,col_both,clusterAnaResults_TrainVsGen)

%  Session comparison- Cluster analyses
for bComp=1:2
    if bComp==1
        col=[c_t(1,:);c_t(2,:)];
        sha_clamp_1_2={ha.TC_tw_C1; ha.TC_tw_C2};
        clusterAnaResults_Clamp=clusterAnaResults_Clamp_t1_t2;
    elseif bComp==2
        col=[c_t(1,:);c_g(2,:)];
        sha_clamp_1_2={ha.TC_tw_C1; ha.TC_gw_C2};
        clusterAnaResults_Clamp=clusterAnaResults_Clamp_t1_g2;
    end

    plotMeanTimeCourse_Interference_SepWedges_clusterAna(sha_clamp_1_2,prot,nC_base,clusterAnaResults_Clamp,col)
end


%% Summary analysis plots

% aftereffect - all cycles
yLabel='Reach angle (deg)';
yLim=[5 20; -6 6; -24 24];
yTick={0:5:50; -6:6:6; -20:10:20};
col=[c_t(1,:);c_t(2,:);c_g(2,:)];
plotSumAnalysis_3within_bar_violin_indiv_diff_AntiCLongV(ha.mAftereffect_all,col,yLabel,yLim,yTick)

