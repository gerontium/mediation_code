% 04/02/15 THIS IS THE ONE
clear all
close all
clc
pause(0.5)

save_path = 'C:\Users\Mick\Documents\Ger\Evidence Acculumation Project\Writing\Journal Article\Plots for paper\05_11_14\';

disp('Loading...')
load('jim_corr_ref2','all_RT_time8Hz_200ms_20ms_Pz_POz','allsubj')
load('corrs_35Hz_filt_200ms','all_stim_time_CPP','all_resp_time_CPP','all_stim_time_peak','smooth_time_stim_amp','smooth_time_stim_slope', ...
    'smooth_time_resp_amp','smooth_time_resp_slope')

disp('Loaded')

all_RT_time = all_RT_time8Hz_200ms_20ms_Pz_POz;
all_stim_time = all_stim_time_CPP;
all_resp_time = all_resp_time_CPP;

exclude_chans = [24 28 33 37 47 61];
plot_chans = [1:23,25:27,29:32,34:36,38:46,48:60,62:64];
% exclude_chans = [];
% plot_chans = [1:64];
left_hemi = [1:23,25:27];
right_hemi = [34:36,39:46,49:60,62:64];
centre_chans = [29:32,38,48];
posterior_chans = [16:32,53:64]; posterior_chans(find(ismember(posterior_chans,exclude_chans))) = [];
alpha_chans = [25,27,62,64];

ch_lr = [60;23];
% ch_rl = [16;53];
ch_rl = [23;60];
ch_CPP = [31];
timer = [256,356,246];

plot_titles = {'Contra N2pc','Ipsi N2pc','CMI N2pc','CPP amplitude','CPP slope','LRP amplitude','LRP slope'};
% plot_titles = {'Contra N2pc','Ipsi N2pc','CMI N2pc'};

yellow = [250,250,25]; blue = [25,50,250]; purple = [250,25,250]; red = [250,25,25];
green = [50,250,25]; orange = [250,135,25]; cyan = [25,250,250];

nice_colours = {blue,red,green,cyan,purple,orange,yellow};

counter = 1;
for s = 1:length(allsubj)
    all_RT_time_group(1,counter:counter+size(all_RT_time{s},2)-1) = single(all_RT_time{s});
    all_stim_time_group(:,:,counter:counter+size(all_stim_time{s},3)-1) = single(all_stim_time{s});
    all_resp_time_group(:,:,counter:counter+size(all_resp_time{s},3)-1) = single(all_resp_time{s});
    all_stim_time_peak_group(:,counter:counter+size(all_stim_time_peak{s},2)-1) = single(all_stim_time_peak{s});
    counter = counter+size(all_RT_time{s},2);
end

%% mediation for posh plots
% x = N2
% y = RT
% m = CPP
% 1 a   X -> M relationship
% 2 b   M -> Y relationship
% 3 cp  unmediated X -> Y relationship (residual)
% 4 c   X -> Y relationship
% 5 ab  mediated X -> Y by M (a * b)
clear stim_ab
for tt = 1:size(all_stim_time_group,2)
    disp(tt)
    [paths, toplevelstats, firstlevelstats] = mediation(squeeze(all_stim_time_peak_group(1,:))', ...
        squeeze(all_RT_time_group(1,:))',squeeze(all_stim_time_group(1,tt,:)),'stats','covs',squeeze(all_stim_time_peak_group(3,:))');
    stim_ab(tt,:,1) = toplevelstats.mean;
    stim_ab(tt,:,2) = toplevelstats.p;
end
figure
subplot(2,1,1)
plot(smooth_time_stim_slope,squeeze(stim_ab(:,5,1)),'LineWidth',2), hold on
line([xlim],[0,0],'Color','k','LineWidth',1.5);
title('Mediation effect across time of stim-locked CPP slope on N2c-RT relationship','FontSize',15)
set(gca,'FontSize',20);
subplot(2,1,2)
plot(smooth_time_stim_slope,squeeze(stim_ab(:,5,2)),'LineWidth',2), hold on
line([xlim],[0.05,0.05],'Color','k','LineWidth',1.5,'LineStyle','--');
set(gca,'FontSize',20);
% return
% close all
% for i = 1:5
%     figure
%     plot(smooth_time_stim_slope,squeeze(stim_ab(:,i,1)))
% end
% figure
% plot(smooth_time_stim_slope,squeeze(stim_medi(:,5,2)))

conder = 1;
clear extreme_ab
N2pc_vec_orig = squeeze(all_stim_time_peak_group(conder,:))';
CPP_vec_orig = squeeze(all_stim_time_peak_group(3,:))';
RT_vec_orig = squeeze(all_RT_time_group(1,:))';
for i = 1:1000
    disp(i)
    rander = randi(length(N2pc_vec_orig),[1,length(N2pc_vec_orig)]);
    N2pc_vec = N2pc_vec_orig(rander);
    clear ab_rand
    for tt = 1:find(smooth_time_stim_slope==601)%size(all_resp_time_group,2)
        CPP_vec = squeeze(all_stim_time_group(1,tt,:));
        [~, toplevelstats, ~] = mediation(N2pc_vec, ...
            RT_vec_orig,CPP_vec,'stats','covs',CPP_vec_orig);
        ab_rand(tt,1) = toplevelstats.mean(5);
    end
    max_ab_temp = max(ab_rand);
    min_ab_temp = min(ab_rand);
    if max_ab_temp>abs(min_ab_temp)
        extreme_ab(i) = max_ab_temp;
    else
        extreme_ab(i) = min_ab_temp;
    end
    clc
end

figure
hist(extreme_ab);

% prc_temp = prctile(extreme_beta,[5 95])
prc_temp_stim = prctile(extreme_ab,[2.5 97.5])
prc_temp_stim2 = prctile(extreme_ab,[1 99])
prc_temp_stim3 = prctile(extreme_ab,[0.1 99.9])

coloro = {'b','r','c'}; timero = [0.04,0.05,0.06];
clear h h2
figure
h(1) = plot(smooth_time_stim_slope,squeeze(stim_ab(:,5,1)),'b','LineWidth',5); hold on

set(gca,'FontSize',24,'xlim',[100,600],'xtick',[100:100:600]);%,'ylim',[-0.05,0.05],'ytick',[-0.04:0.02:0.04]);
xlabel('Time (ms)','FontName','Arial','FontSize',24)
ylabel('Correlation with N2pc amplitude (Pearson''s R)','FontName','Arial','FontSize',20)
line([xlim],[0,0],'Color','k','LineWidth', 1.4);
line([xlim],[prc_temp_stim(1),prc_temp_stim(1)],'Color','k','LineWidth', 2,'LineStyle','--');
line([xlim],[prc_temp_stim(2),prc_temp_stim(2)],'Color','k','LineWidth', 2,'LineStyle','--');
line([xlim],[prc_temp_stim2(1),prc_temp_stim2(1)],'Color','k','LineWidth', 2,'LineStyle','--');
line([xlim],[prc_temp_stim2(2),prc_temp_stim2(2)],'Color','k','LineWidth', 2,'LineStyle','--');
line([xlim],[prc_temp_stim3(1),prc_temp_stim3(1)],'Color','k','LineWidth', 2,'LineStyle','--');
line([xlim],[prc_temp_stim3(2),prc_temp_stim3(2)],'Color','k','LineWidth', 2,'LineStyle','--');
legend(h,{'Contra N2pc'}, ...
    'Location','SouthEast');
save([save_path,'stimlocked_N2c_CPP_med.mat'],'smooth_time_stim_slope','stim_ab','prc_temp_stim','extreme_ab')
disp('done')
% return
%% mediation for posh plots - response locked
% x = N2
% y = RT
% m = CPP
% 1 a   X -> M relationship
% 2 b   M -> Y relationship
% 3 cp  unmediated X -> Y relationship (residual)
% 4 c   X -> Y relationship
% 5 ab  mediated X -> Y by M (a * b)
clear resp_ab
for tt = 1:size(all_resp_time_group,2)
    disp(tt)
    [paths, toplevelstats, firstlevelstats] = mediation(squeeze(all_stim_time_peak_group(1,:))', ...
        squeeze(all_RT_time_group(1,:))',squeeze(all_resp_time_group(1,tt,:)),'stats','covs',squeeze(all_stim_time_peak_group(3,:))');
    resp_ab(tt,:,1) = toplevelstats.mean;
    resp_ab(tt,:,2) = toplevelstats.p;
end
figure
subplot(2,1,1)
plot(smooth_time_resp_slope,squeeze(resp_ab(:,5,1)),'LineWidth',2), hold on
line([xlim],[0,0],'Color','k','LineWidth',1.5);
title('Mediation effect across time of resp-locked CPP slope on N2c-RT relationship','FontSize',15)
set(gca,'FontSize',20);
subplot(2,1,2)
plot(smooth_time_resp_slope,squeeze(resp_ab(:,5,2)),'LineWidth',2), hold on
line([xlim],[0.05,0.05],'Color','k','LineWidth',1.5,'LineStyle','--');
set(gca,'FontSize',20);

% close all
% for i = 1:5
%     figure
%     plot(smooth_time_resp_slope,squeeze(resp_medi(:,i,1)))
% end
% figure
% plot(smooth_time_resp_slope,squeeze(resp_medi(:,5,2)))

conder = 1;
clear extreme_ab
N2pc_vec_orig = squeeze(all_stim_time_peak_group(conder,:))';
CPP_vec_orig = squeeze(all_stim_time_peak_group(3,:))';
RT_vec_orig = squeeze(all_RT_time_group(1,:))';
for i = 1:1000
    disp(i)
    rander = randi(length(N2pc_vec_orig),[1,length(N2pc_vec_orig)]);
    N2pc_vec = N2pc_vec_orig(rander);
    clear ab_rand
    for tt = find(smooth_time_resp_slope==-501):find(smooth_time_resp_slope==-101)
        CPP_vec = squeeze(all_resp_time_group(1,tt,:));
        [~, toplevelstats, ~] = mediation(N2pc_vec, ...
            RT_vec_orig,CPP_vec,'stats','covs',CPP_vec_orig);
        ab_rand(tt,1) = toplevelstats.mean(5);
    end
    max_ab_temp = max(ab_rand);
    min_ab_temp = min(ab_rand);
    if max_ab_temp>abs(min_ab_temp)
        extreme_ab(i) = max_ab_temp;
    else
        extreme_ab(i) = min_ab_temp;
    end
    clc
end

figure
hist(extreme_ab);

% prc_temp = prctile(extreme_beta,[5 95])
prc_temp_resp = prctile(extreme_ab,[2.5 97.5])
prc_temp_resp2 = prctile(extreme_ab,[1 99])
prc_temp_resp3 = prctile(extreme_ab,[0.1 99.9])

coloro = {'b','r','c'}; timero = [0.04,0.05,0.06];
clear h h2
figure
h(1) = plot(smooth_time_resp_slope,squeeze(resp_ab(:,5,1)),'b','LineWidth',5); hold on

set(gca,'FontSize',24,'xlim',[-600,0]);%,'ylim',[-0.05,0.05],'ytick',[-0.04:0.02:0.04]);
xlabel('Time (ms)','FontName','Arial','FontSize',24)
ylabel('Correlation with N2pc amplitude (Pearson''s R)','FontName','Arial','FontSize',20)
line([xlim],[0,0],'Color','k','LineWidth', 1.4);
% line([xlim],[prc_temp_resp(1),prc_temp_resp(1)],'Color','k','LineWidth', 2,'LineStyle','--');
% line([xlim],[prc_temp_resp(2),prc_temp_resp(2)],'Color','k','LineWidth', 2,'LineStyle','--');
% line([xlim],[prc_temp_resp2(1),prc_temp_resp2(1)],'Color','k','LineWidth', 2,'LineStyle','--');
% line([xlim],[prc_temp_resp2(2),prc_temp_resp2(2)],'Color','k','LineWidth', 2,'LineStyle','--');
% line([xlim],[prc_temp_resp3(1),prc_temp_resp3(1)],'Color','k','LineWidth', 2,'LineStyle','--');
% line([xlim],[prc_temp_resp3(2),prc_temp_resp3(2)],'Color','k','LineWidth', 2,'LineStyle','--');
legend(h,{'Contra N2pc'}, ...
    'Location','SouthEast');
save([save_path,'resplocked_N2c_CPP_med.mat'],'smooth_time_resp_slope','resp_ab','prc_temp_resp','extreme_ab')




%%
%  ***********************************************************************
%%
clear all
close all
clc
disp('Loading...')
load('centre_corr_peak')
disp('Loaded')

save_path = 'C:\Users\Mick\Documents\Ger\Evidence Acculumation Project\Writing\Journal Article\Plots for paper\05_11_14\';

% ************************ GENERAL PLOTS VS RT ****************************
all_RT_time{s}(1,c_total_t:c_total_t+length(RT_temp)-1) = zscore(log(RT_temp));

all_stim_time_CPP{s}(1,:,c_total_t:c_total_t+size(CPP_stim_slope,2)-1) = zscore(CPP_stim_slope,0,2);
all_resp_time_CPP{s}(1,:,c_total_t:c_total_t+size(CPP_resp_slope,2)-1) = zscore(CPP_resp_slope,0,2);

all_stim_time_peak{s}(1,c_total_t:c_total_t+length(both_hemi_onset_stim_amp)-1) = zscore(both_hemi_onset_stim_amp);
all_stim_time_peak{s}(2,c_total_t:c_total_t+length(CPP_stim_onset_peak_amp)-1) = zscore(CPP_stim_onset_peak_amp);
            
counter = 1;
for s = 1:length(allsubj)
    disp(s)
    all_RT_time_group(1,counter:counter+size(all_RT_time{s}(1,:),2)-1) = all_RT_time{s}(1,:); % RT
    
    all_stim_time_peak_group(:,counter:counter+size(all_stim_time_peak{s},2)-1) = all_stim_time_peak{s}; % left and right hemi N2, CPP also
    
    all_stim_time_CPP_group(1,:,counter:counter+size(all_stim_time_CPP{s},3)-1) = all_stim_time_CPP{s}; % stim CPP slope
    all_resp_time_CPP_group(1,:,counter:counter+size(all_resp_time_CPP{s},3)-1) = all_resp_time_CPP{s}; % resp CPP slope
    counter = counter+size(all_RT_time{s}(1,:),2);
end

plot_titles = {'Right & Left Hemi N2p Onset'};

%% mediation for posh plots
% x = N2
% y = RT
% m = CPP
% 1 a   X -> M relationship
% 2 b   M -> Y relationship
% 3 cp  unmediated X -> Y relationship (residual)
% 4 c   X -> Y relationship
% 5 ab  mediated X -> Y by M (a * b)

clear stim_ab
for tt = 1:size(all_stim_time_CPP_group,2)
    disp(tt)
    [paths, toplevelstats, firstlevelstats] = mediation(squeeze(all_stim_time_peak_group(1,:))', ...
        squeeze(all_RT_time_group(1,:))',squeeze(all_stim_time_CPP_group(1,tt,:)),'stats','covs',squeeze(all_stim_time_peak_group(2,:))');
    stim_ab(tt,:,1) = toplevelstats.mean;
    stim_ab(tt,:,2) = toplevelstats.p;
end
figure
subplot(2,1,1)
plot(smooth_time_stim_slope,squeeze(stim_ab(:,4,1)),'LineWidth',2)
subplot(2,1,2)
plot(smooth_time_stim_slope,squeeze(stim_ab(:,4,2)),'LineWidth',2)

line([xlim],[0,0],'Color','k','LineWidth',1.5);
title('Mediation effect across time of stim-locked CPP slope on N2b-RT relationship','FontSize',15)
set(gca,'FontSize',20);
subplot(2,1,2)
plot(smooth_time_stim_slope,squeeze(stim_ab(:,5,2)),'LineWidth',2), hold on
line([xlim],[0.05,0.05],'Color','k','LineWidth',1.5,'LineStyle','--');
set(gca,'FontSize',20);

conder = 1;
clear extreme_ab
N2pc_vec_orig = squeeze(all_stim_time_peak_group(conder,:))';
CPP_vec_orig = squeeze(all_stim_time_peak_group(2,:))';
RT_vec_orig = squeeze(all_RT_time_group(1,:))';
for i = 1:1000
    disp(i)
    rander = randi(length(N2pc_vec_orig),[1,length(N2pc_vec_orig)]);
    N2pc_vec = N2pc_vec_orig(rander);
    clear ab_rand
    for tt = find(smooth_time_stim_slope<=700)%size(all_resp_time_group,2)
        CPP_vec = squeeze(all_stim_time_CPP_group(1,tt,:));
        [~, toplevelstats, ~] = mediation(N2pc_vec, ...
            RT_vec_orig,CPP_vec,'stats','covs',CPP_vec_orig);
        ab_rand(tt,1) = toplevelstats.mean(5);
    end
    max_ab_temp = max(ab_rand);
    min_ab_temp = min(ab_rand);
    if max_ab_temp>abs(min_ab_temp)
        extreme_ab(i) = max_ab_temp;
    else
        extreme_ab(i) = min_ab_temp;
    end
    clc
end

figure
hist(extreme_ab);

% prc_temp = prctile(extreme_beta,[5 95])
prc_temp_stim = prctile(extreme_ab,[2.5 97.5])
prc_temp_stim2 = prctile(extreme_ab,[1 99])
prc_temp_stim3 = prctile(extreme_ab,[0.1 99.9])

clear h h2
figure 
h(1) = plot(smooth_time_stim_slope,squeeze(stim_ab(:,5,1)),'b','LineWidth',5); hold on
set(gca,'FontSize',20,'xlim',[100,600],'ylim',[-0.04,0.04]);
xlabel('Time (ms)','FontName','Arial','FontSize',20)
ylabel('Indirect effect ''ab'' of N2b on RT via CPP slope','FontName','Arial','FontSize',20)

line([xlim],[0,0],'Color','k','LineWidth', 1.5);

h2(1) = line([xlim],[prc_temp_stim(1),prc_temp_stim(1)],'Color','k','LineWidth', 2,'LineStyle','--');
line([xlim],[prc_temp_stim(2),prc_temp_stim(2)],'Color','k','LineWidth', 2,'LineStyle','--');
legend(h,{'Left & Right Hemi N2b'}, ...
    'Location','NorthWest');
save([save_path,'stimlocked_N2b_CPP_med.mat'],'smooth_time_stim_slope','stim_ab','prc_temp_stim','extreme_ab')
% % load([save_path,'stimlocked_N2b_CPP_med.mat'],'smooth_time_stim_slope','stim_ab','prc_temp_stim','extreme_ab')
%% mediation for posh plots - response locked
% x = N2
% y = RT
% m = CPP
% 1 a   X -> M relationship
% 2 b   M -> Y relationship
% 3 cp  unmediated X -> Y relationship (residual)
% 4 c   X -> Y relationship
% 5 ab  mediated X -> Y by M (a * b)

clear resp_ab
for tt = 1:size(all_resp_time_CPP_group,2)
    disp(tt)
    [paths, toplevelstats, firstlevelstats] = mediation(squeeze(all_stim_time_peak_group(1,:))', ...
        squeeze(all_RT_time_group(1,:))',squeeze(all_resp_time_CPP_group(1,tt,:)),'stats','covs',squeeze(all_stim_time_peak_group(2,:))');
    resp_ab(tt,:,1) = toplevelstats.mean;
    resp_ab(tt,:,2) = toplevelstats.p;
end
figure
subplot(2,1,1)
plot(smooth_time_resp_slope,squeeze(resp_ab(:,5,1)),'LineWidth',2), hold on
line([xlim],[0,0],'Color','k','LineWidth',1.5);
title('Mediation effect across time of resp-locked CPP slope on N2b-RT relationship','FontSize',15)
set(gca,'FontSize',20);
subplot(2,1,2)
plot(smooth_time_resp_slope,squeeze(resp_ab(:,5,2)),'LineWidth',2), hold on
line([xlim],[0.05,0.05],'Color','k','LineWidth',1.5,'LineStyle','--');
set(gca,'FontSize',20);

conder = 1;
clear extreme_ab
N2pc_vec_orig = squeeze(all_stim_time_peak_group(conder,:))';
CPP_vec_orig = squeeze(all_stim_time_peak_group(2,:))';
RT_vec_orig = squeeze(all_RT_time_group(1,:))';
for i = 1:1000
    disp(i)
    rander = randi(length(N2pc_vec_orig),[1,length(N2pc_vec_orig)]);
    N2pc_vec = N2pc_vec_orig(rander);
    clear ab_rand
    for tt = length(find(smooth_time_resp_slope<=-600)):length(find(smooth_time_resp_slope<=-100))%size(all_resp_time_group,2)
        CPP_vec = squeeze(all_resp_time_CPP_group(1,tt,:));
        [~, toplevelstats, ~] = mediation(N2pc_vec, ...
            RT_vec_orig,CPP_vec,'stats','covs',CPP_vec_orig);
        ab_rand(tt,1) = toplevelstats.mean(5);
    end
    max_ab_temp = max(ab_rand);
    min_ab_temp = min(ab_rand);
    if max_ab_temp>abs(min_ab_temp)
        extreme_ab(i) = max_ab_temp;
    else
        extreme_ab(i) = min_ab_temp;
    end
    clc
end

figure
hist(extreme_ab);

% prc_temp = prctile(extreme_beta,[5 95])
prc_temp_resp = prctile(extreme_ab,[2.5 97.5])
prc_temp_resp2 = prctile(extreme_ab,[1 99])
prc_temp_resp3 = prctile(extreme_ab,[0.1 99.9])

clear h h2
figure 
h(1) = plot(smooth_time_resp_slope,squeeze(resp_ab(:,5,1)),'b','LineWidth',5); hold on
set(gca,'FontSize',20,'xlim',[-600,0],'ylim',[-0.015,0.015]);
ylabel('Indirect effect ''ab'' of N2b on RT via CPP slope','FontName','Arial','FontSize',20)
xlabel('Time (ms)','FontName','Arial','FontSize',20)
line([xlim],[0,0],'Color','k','LineWidth', 1.5);

h2(1) = line([xlim],[prc_temp_resp(1),prc_temp_resp(1)],'Color','k','LineWidth', 2,'LineStyle','--');
line([xlim],[prc_temp_resp(2),prc_temp_resp(2)],'Color','k','LineWidth', 2,'LineStyle','--');
legend(h,{'Left & Right Hemi N2b'}, ...
    'Location','NorthWest');
save([save_path,'resplocked_N2b_CPP_med.mat'],'smooth_time_resp_slope','resp_ab','prc_temp_resp','extreme_ab')
% % prc_temp_resp = prc_temp_stim;
% % save([save_path,'resplocked_N2b_CPP_med.mat'],'prc_temp_resp','extreme_ab','-append')
return


