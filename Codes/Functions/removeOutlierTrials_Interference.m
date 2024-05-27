function [ha, nT_org, nT_removed] = removeOutlierTrials_Interference(ha_org)

nT_org=length(ha_org);
ha_withOutlierTrials=ha_org;

spikeThre=20;

trial_out_abs=find(abs(ha_withOutlierTrials)>100); % base on absolute hand angle
if abs(ha_withOutlierTrials(1))>spikeThre % when the first trial has too large error
    trial_out_abs=[1 trial_out_abs];
end
ha_withOutlierTrials(trial_out_abs)=nan;

% Spike of 1 trial: a "spike" is when there is an oulier trial in which error was suddenly very different from the trial before and the one / two trials after
dha_s=diff(ha_withOutlierTrials); % base on sudden change in hand angle
diffLarge=find(abs(dha_s)>spikeThre)+1; % find how errors change from trial to trial
iSpike1=find(diff(diffLarge)==1); % 
trial_out1=diffLarge(iSpike1);
ha_withOutlierTrials(trial_out1)=nan;

trial_out=[trial_out_abs trial_out1];
nT_removed=length(trial_out);

% figSize=[50 100 1600 400];
% figure('position',figSize)
% plot(1:nT_org,ha_org,trial_out,ha_org(trial_out),'x')
% pause

ha=ha_withOutlierTrials;


end

