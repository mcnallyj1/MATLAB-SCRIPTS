clear

%% Script for analysis of Go/No-Go LFP and Respiratory Data 
% Updated 9/11/24 JMM
% Allows:
% 1)Spectral analysis (wavelet transform) of total or evoked activity
%   during trial performance
% 2)PAC analysis of LFP data
% 3)Respiratory analysis for mice that have thermistor data

%% Data Entry

% enter excel file name  
filelist = 'Therm22multifile';
mouseID = 'Therm22';

%% Behavioral response of interest

graph = "FALSE ALARMS"; % HITS, MISSES, CORRECT REJECTS, FALSE ALARMS

%% LFP channel to analyze

LFP_Channel = "BF"; % CA1, V1, BF, AI, OB, Prl, RS

%% Analysis to perform?
spectral = 'Y'; % Y or N
respiration = 'N'; % Y or N
PAC = 'N'; % Y or N

t0 = "np_exit"; % np_entry, np_exit, cue_off %% What analysis is time-locked to (denotes time zero)?

%% Type of spectral analysis??

spec_type = 'total'; % evoked or total


%%% time normalization?
time_warp = 'Y'; %% this analysis is only good with np_enter or np_exit time locking
start_scale = 4; % define time before to include
end_scale = 3; % define time after completion to include


%% Filter & remove trial records with artifacts

notch_filter = 'Y'; % Y or N
remove_artifacts = 'Y'; % Y or N
show_outliers = 'N'; % Y or N
percent_sat = 5; % percent of total trial length at or above max physiological voltage
%max_v = 350; % max physiological voltage for channel in uV

%% Spectral Params & Bin Size
pre = 10; % secs
post = 5; % secs
pre_resp = 10; % secs prior to trial end for resp analysis
pre_pac = 10;  % secs prior to trial end for PAC analysis 

smf = 50; %smoothing factor # of pts
freq = [1 100];
dfreq = .2;

%% Downsample data?
down_sample = 'Y';
ds = 4; % downsample data (Fs/ds)

%% Plotting parameters
cmap = 'jet'; %color map for spectral plot 
c_lo = 5;
c_hi = 40;

%% Wavelet parameters
fb = 28; % 13.5 bandwidth
fc = 0.5; % 0.5 center freq

%% Define frequency bands of interest

theta = [6 10];
alpha = [10 15];
beta = [15 30];
gamma = [30 90];
epsilon = [90 200];

%% PAC parameters_
phase_freq = [1 15];
amp_freq = [15 100];

%% import data 

[num, txt] = xlsread([filelist,'.xlsx']);
filesLFP = txt(:,1);
filesTTL = txt(:,2);

rewDisp_channel = 3; % syringe pump 
nosepk_channel = 5;
lick_channel = 6; % head entry to reward port 
cue_channel = 7; % cue light
odor1_channel = 8; % (+) Go stimulus
odor2_channel = 9; % (-) NoGo stimulus
resp_channel = 10; 

data = [];
nosepk_times = [];
nosepk_level = [];
odor1_times = [];
odor2_times = [];
reward_times = [];
retrieval_times = [];
resp_data = [];

endtime = 0;

for y = 1:length(filesLFP)
    fileLFP = load([filesLFP{y}]);
    fileTTL = load([filesTTL{y}]);
    if isfield(fileLFP, 'info_lfp') == 1
        idx = find(strcmp([fileLFP.info_lfp.header.LFP_channel_label], LFP_Channel)); % single line engine
    else
        idx = find(strcmp([fileLFP.info.header.LFP_channel_label], LFP_Channel)); % single line engine
    end
    eval(['pdata = fileLFP.LFP(idx,:);']); % EEG record (mV)
    %eval(['pdata = fileLFP.LFP;']); % EEG record (mV) % JADE
    data = cat(2,data,pdata);
    %eval(['fs = fileLFP.info.header.sampleRate;']); % EEG sampling frequency Hz (2000)
    fs = 2000;
    eval(['pnosepk_times = fileTTL.' filesTTL{y} '_Ch' num2str(nosepk_channel) '.times;']); % time (secs)
    nosepk_times=cat(1,nosepk_times,(pnosepk_times + endtime));
    eval(['pnosepk_level = fileTTL.' filesTTL{y} '_Ch' num2str(nosepk_channel) '.level;']); % binary
    nosepk_level=cat(1,nosepk_level,pnosepk_level);
    eval(['podor1_times = fileTTL.' filesTTL{y} '_Ch' num2str(odor1_channel) '.times;']); % time (secs)
    odor1_times=cat(1,odor1_times,(podor1_times + endtime));
    eval(['podor2_times = fileTTL.' filesTTL{y} '_Ch' num2str(odor2_channel) '.times;']); % time (secs)
    odor2_times=cat(1,odor2_times,(podor2_times + endtime));
    eval(['preward_times = fileTTL.' filesTTL{y} '_Ch' num2str(rewDisp_channel) '.times;']); % time (secs)
    reward_times=cat(1,reward_times,(preward_times + endtime));
    eval(['pretrieval_times = fileTTL.' filesTTL{y} '_Ch' num2str(lick_channel) '.times;']); % time (secs)
    retrieval_times=cat(1,retrieval_times,(pretrieval_times + endtime));
    eval(['presp_data = fileTTL.' filesTTL{y} '_Ch' num2str(resp_channel) '.values*1000;']); % time (secs)
    resp_data = cat(1,resp_data,presp_data);
    endtime = length(data)*(1/fs);

    % data = detrend(data);
    % 
    % [A,S] = jade(data);
    % idx = find(strcmp([fileLFP.info_lfp.header.LFP_channel_label], LFP_Channel)); % Finds labels for each LFP channel
    % data = S(idx,:); % Focuses data on channel of interest
    % 
    % tapers = [5 9]; % # of tapers and bandwith = (2*tapers) -1)
    % fpass = [0 100]; %frequency band of interest
    % moving_win = [30 1]; % moving window for spectral analysis [window size, step size] in seconds
    % db = 'l'; %scale plot in db 'l', linear 'n'
    % 
    % params.Fs = fs;
    % params.fpass = fpass;
    % params.tapers = tapers;
    % params.trialave = 1;
    % params.err = 0;
    % 
    % [spec,T,F] = mtspecgramc(data, moving_win, params);
    % 
    % figure;  
    %          plot_matrix(spec,T, F, db);
    %          eval(['colormap ' cmap ';']);
    %          %caxis([-40 -15]);
    %          xlabel('Time (s)');
    %          ylim([0 100]);
    %          ylabel('Frequency (Hz)');
    %          title('0 mg/kg', LFP_Channel);

%  [S,f]=mtspectrumc(data,params);
% plot(f,pow2db(S))
% grid on
% xlabel('Frequency (Hz)')
% ylabel('Power Spectrum (dB)')
% title('Power Spectra')
end


%% Clean up data

data = detrend(data);
max_v = 3*rms(data);

% 60/120/180 Hz notch filters
if notch_filter == 'Y'
    d1 = designfilt('bandstopiir','FilterOrder',10, ...
               'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5, ...
               'DesignMethod','butter','SampleRate',fs);

     d2 = designfilt('bandstopiir','FilterOrder',10, ...
                'HalfPowerFrequency1',119.5,'HalfPowerFrequency2',120.5, ...
                'DesignMethod','butter','SampleRate',fs);
     
     d3 = designfilt('bandstopiir','FilterOrder',10, ...
                'HalfPowerFrequency1',179.5,'HalfPowerFrequency2',180.5, ...
                'DesignMethod','butter','SampleRate',fs);

    data = filtfilt(d1,data);
    % data = filtfilt(d2,data);
    data = filtfilt(d3,data);
end

% downsample data and time
if down_sample == 'Y'
    data = resample(data,fs/ds,fs);
    resp_data = resample(resp_data,fs/ds,fs);
    fs = fs/ds; % convert fs to match 
end


%% Derive scales for wavelet analysis
freq2 = freq(1):dfreq:freq(2);
scales = frq2scale(freq2,['cmor' num2str(fb) '-' num2str(fc)], 1/fs);

%% Make table of trial times (nose poke times) for all trials 

% calculate nose poke duration based on odor times
for i = 1:2:length(odor2_times)-1
    NOGOnosePokeDur(i) = odor2_times(i+1) - odor2_times(i);
end

for i = 1:2:length(odor1_times)-1
    GOnosePokeDur(i) = odor1_times(i+1) - odor1_times(i);
end

% select nose pokes that are greater than 200ms
% then find the nosePoke data points that correspond with odor delivery

NOGOnosePokes = find(NOGOnosePokeDur > 0.200); 
for i = 1:length(NOGOnosePokes)
    x = find(nosepk_times < odor2_times(NOGOnosePokes(i)) & nosepk_level == 1, 1, "last");
    NOGOlabel(i,1) = "NOGO";
    NOGOtrials(i,1) = nosepk_times(x);
    NOGOtrials(i,2) = nosepk_times(x+1);
end

GOnosePokes = find(GOnosePokeDur > 0.200); 
for i = 1:length(GOnosePokes)
    x = find(nosepk_times < odor1_times(GOnosePokes(i)) & nosepk_level == 1, 1, "last");
    GOlabel(i,1) = "GO";
    GOtrials(i,1) = nosepk_times(x);
    GOtrials(i,2) = nosepk_times(x+1);
end

NOGO = table(NOGOlabel, NOGOtrials(:,1), NOGOtrials(:,2), 'VariableNames', {'TrialType', 'start','end'});
GO = table(GOlabel, GOtrials(:,1), GOtrials(:,2), 'VariableNames', {'TrialType', 'start','end'});
trials = [GO; NOGO];
trials = sortrows(trials, "start");
premature = (height(odor1_times) + height(odor2_times)) - height(trials); % just a number for now, can find where premature trials are later

%% Determine and label behavioral responses for each trial in table

for i = 1:size(trials)
   if trials.TrialType(i) == "GO"
       x = find(reward_times > trials.end(i), 1, "first");

       % if reward is dispensed or not
       if reward_times(x) < trials.end(i)+3
           trials.Response(i) = "hit";
       else
           trials.Response(i) = "miss";
       end

   elseif trials.TrialType(i) == "NOGO"
       x = find(retrieval_times > trials.end(i), 1, "first");

       % if trial is re-initiated or if the mouse does nothing for 3sec, otherwise incorrect
       if (i ~= size(trials,1) && trials.start(i+1) < trials.end(i)+3) || (retrieval_times(x) > trials.end(i)+3)
           trials.Response(i) = "correctReject";
       else
           trials.Response(i) = "falseAlarm";
       end
   end
end

%% determine amount of time between trials

trials.delta(1) = trials.start(1);

for x = 2:size(trials,1)
    trials.delta(x) = trials.start(x)-trials.end(x-1);
end

%% Check number of trials

idx_hit = find(trials.Response == "hit");
idx_miss = find(trials.Response == "miss");
idx_corrRej = find(trials.Response == "correctReject");
idx_falseAlarm = find(trials.Response == "falseAlarm");

%%% record trial length
for x = 1:size(trials)
    trials.length(x) = trials.end(x)-trials.start(x);
end

%% Determine the Reaction time for each trial in table (time from odor cue onset to response)

RT_hit = zeros(size(idx_hit,1),1);

for x = 1: size(idx_hit,1)
    y = find(reward_times > trials.end(idx_hit(x)), 1, "first");
    RT_hit(x) = reward_times(y)- trials.start(idx_hit(x)); 
    z = find(trials.end > trials.start(idx_hit(x)), 1, "first");
    RT2_hit(x) = trials.end(z) - trials.start(idx_hit(x));
end

for x = 1: size(idx_corrRej,1)
    y = find(nosepk_times > trials.end(idx_corrRej(x)), 1, "first");
    RT_CR(x) = nosepk_times(y)- trials.start(idx_corrRej(x)); 
    z = find(trials.end > trials.start(idx_corrRej(x)), 1, "first");
    RT2_CR(x) = trials.end(z) - trials.start(idx_corrRej(x));
end

for x = 1: size(idx_falseAlarm,1)
    y = find(reward_times > trials.end(idx_falseAlarm(x)), 1, "first");
    RT_FA(x) = reward_times(y)- trials.start(idx_falseAlarm(x)); 
    z = find(trials.end > trials.start(idx_falseAlarm(x)), 1, "first");
    RT2_FA(x) = trials.end(z) - trials.start(idx_falseAlarm(x));
end

aveRT_hit = mean(RT_hit,1);
aveRT2_hit = mean(RT2_hit,1);

aveRT_CR = mean(RT_CR,1);
aveRT2_CR = mean(RT2_CR,1);

aveRT_FA = mean(RT_FA,1);
aveRT2_FA = mean(RT2_FA,1);



%% Extract EEG 

% For each behavioral type (hit, miss, falseAlarm, correctReject)
switch t0
    case "np_entry"
        type = 'trials.start';            
    case "np_exit"
        type = 'trials.end';   
    case "cue_off" 
        type = 'cue.end';        
end

hits = [];
misses = [];
correctRejects = [];
falseAlarms = [];

% Assign trial_LFP to behavioral response of interest
h = 1;
m = 1;
f = 1;
c = 1;

for i = 1:size(trials,1); % Add - 1 outside of the parentheses if that one error happens
    if trials.Response(i) == "hit"
        eval(['hits(:,h) = data((' type '(i) - pre)*fs : (' type '(i) + post)*fs);']) % Table where columns = trials
        hit_length(h) = trials.length(i);
        h = h + 1;
    elseif trials.Response(i) == "correctReject"
        eval(['correctRejects(:,c) = data((' type '(i) - pre)*fs : (' type '(i) + post)*fs);'])
        CR_length(c) = trials.length(i);
        c = c + 1;
    elseif trials.Response(i) == "miss"
        eval(['misses(:,m) = data((' type '(i) - pre)*fs : (' type '(i) + post)*fs);'])
        miss_length(m) = trials.length(i);
        m = m + 1;
    elseif trials.Response(i) == "falseAlarm"
        eval(['falseAlarms(:,f) = data((' type '(i) - pre)*fs : (' type '(i) + post)*fs);'])
        FA_length(f) = trials.length(i);
        f = f + 1;
    end
end

trial_LFP = [];

switch graph
    case "HITS"
        trial_LFP = hits;
        trial_length = hit_length;
    case "CORRECT REJECTS"
        trial_LFP = correctRejects;
        trial_length = CR_length;
    case "MISSES"
        trial_LFP = misses;
        trial_length = miss_length;
    case "FALSE ALARMS"
        trial_LFP = falseAlarms;
        trial_length = FA_length;
end

%% Artifact filtering

a = 1; % counter for max voltage artifacts

if remove_artifacts == 'Y'

    sf = []; 
    ent = []; 
    max_pts = []; 
    maxx = []; 
    activity = []; 
    mobility =[]; 
    complexity = [];

    for i = 1:size(trial_LFP,2)
        sf(i)= sum(spectralFlatness(trial_LFP(:,i),fs));
        [activity(i), mobility(i), complexity(i)] = hjorth(trial_LFP(:,i));
        ent(i) = mean(wentropy(trial_LFP(:,i)'));
        %FNN(i) = f_fnn(trial_LFP(:,i),2,1,15,2);
        b =  find(abs(trial_LFP(:,i)) >= max_v);
        max_pts(i) = length(b)/size(trial_LFP,1)*100;
        if length(b)/size(trial_LFP,1) > percent_sat/100 % percentage of saturation above which trial will be omitted from spectral analysis
             maxx(a) = i;
             a = a + 1;
        end
    end

    if show_outliers == 'Y'

        figure;
        plot(sf, 'b.')
        title('Spectral Flatness')
        
        figure;
        plot(activity, 'b.')
        title('Hjorth Activity')
        
        figure;
        plot(mobility, 'b.')
        title('Hjorth Mobility')
        
        figure;
        plot(complexity, 'b.')
        title('Hjorth Complexity')
        
        figure;
        plot(ent, 'b.')
        title('Shannon Entropy')

        figure;
        plot(max_pts, 'b.')
        title('Percent of Data Points Above Max Voltage')

    end
    
    %% Suggest outliers

    ave_sf = mean(sf);
    ave_activity = median(activity);
    ave_mobility = median(mobility);
    ave_complexity = median(complexity);
    ave_ent = median(ent);
    
    SD_sf = std(sf);
    SD_activity = std(activity);
    SD_mobility = std(mobility);
    SD_complexity = std(complexity);
    SD_ent = std(ent);
    
    err_sf = find(sf > ave_sf + 2*SD_sf | sf < ave_sf - 2*SD_sf);
    err_activity = find(activity > ave_activity + 2*SD_activity | activity < ave_activity - 2*SD_activity);
    err_mobility = find(mobility > ave_mobility + 2*SD_mobility | mobility < ave_mobility - 2*SD_mobility);
    err_complexity = find(complexity > ave_complexity + 2*SD_complexity | complexity < ave_complexity - 2*SD_complexity);
    err_ent = find(ent > ave_ent + 2*SD_ent | ent < ave_ent - 2*SD_ent);

    if show_outliers == 'Y'
        err_sf
        err_activity
        err_complexity
        err_ent
        maxx
    end
    
    %figure;
    %plot(FNN, 'b.')
    %title('False Nearest Neighbor')
    bad_trials = [err_sf, err_activity, err_mobility, err_complexity, err_ent, maxx];
    artifacts = unique(bad_trials);

    trial_LFP(:,artifacts) = [];
    trial_length(artifacts)= [];
end

%% Plot Total Power

clear pdata

if spectral == 'Y'
switch spec_type
    case 'total'

        if length(trial_LFP) > 1
            for x = 1:size(trial_LFP,2)
                [a_cwt(:,:,x),mor_hz] = cwt(trial_LFP(:,x),scales,['cmor' num2str(fb) '-' num2str(fc)], 1/fs); % Wavelet transform 
                trialsSpec(:,:,x) = abs(a_cwt(:,:,x)).^2;
            end
        end
    
    case 'evoked'
        [a_cwt,mor_hz] = cwt(mean(trial_LFP,2),scales,['cmor' num2str(fb) '-' num2str(fc)], 1/fs); % Wavelet transform 
                trialsSpec = abs(a_cwt).^2;
end

time = (1/fs:1/fs:(pre + post))-pre;

    if time_warp == 'Y'
        num_pts = ((start_scale*100)+(end_scale*100)+100)/0.2;
        dx = linspace(-start_scale*100, end_scale*100+100, num_pts)+start_scale*100;


       for x = 1:size(trial_LFP,2)
            a = trial_length(x); % length from start to end of trial
            for y = 1:size(trialsSpec,1)
                switch t0
                    case "np_entry"
                        Tstart = (pre)*fs-(start_scale*a)*fs;
                        Tend = (pre+a)*fs+(end_scale*a)*fs;
                    case "np_exit"
                        Tstart = (pre-a)*fs-(start_scale*a)*fs;
                        Tend = (pre)*fs+(end_scale*a)*fs;
                end
                interpSpec(y,:,x) = ScaleTime(trialsSpec(y,:,x), Tstart, Tend, num_pts);
            end
        end
        clear trialsSpec

        trialsSpec = interpSpec;

        time = dx-start_scale*100;

    end

    
    TP = mean(trialsSpec,3);

    clear a_cwt 
    
    % 1/f compensation
        for y = 1:size(TP,2)
            TP_fcomp(:,y) = TP(:,y).*mor_hz';
        end

    
 
    % Time-frequency spectrogram
        
    figure;
        imagesc(time, sort(mor_hz,'ascend'), (db(TP_fcomp)));
        eval(['colormap ' cmap ';']);
        set(gca,'YDir','normal');
        hold on
        %xline(0.2, '--','LineWidth',2,'Label','ODOR DELIVERY', 'Color', 'black');
        xline(0, '--','LineWidth',2,'Label','NOSE POKE EXIT', 'LabelHorizontalAlignment', 'Left', 'Color', 'black');
        xlabel('Time (seconds)');
        ylabel('Frequency (Hz)');
        %caxis([65 100])
        title((string(mouseID)) + ' ' + LFP_Channel + ' GNG ' + spec_type + ' Power', graph + ' (' + num2str(size(trial_LFP,2)) + ' trials)');

    % Normalized power plots

    theta_index = find(mor_hz > theta(1) & mor_hz < theta(2));
    alpha_index = find(mor_hz > alpha(1) & mor_hz < alpha(2));
    beta_index = find(mor_hz > beta(1) & mor_hz < beta(2));
    gamma_index = find(mor_hz > gamma(1) & mor_hz < gamma(2));
    epsilon_index = find(mor_hz > epsilon(1) & mor_hz < epsilon(2));
    norm_index = find(time <= pre-2);

    for i = 1:size(time,2)
        ThP(i) = sum(TP_fcomp(theta_index,i));
        AP(i) = sum(TP_fcomp(alpha_index,i));
        BP(i) = sum(TP_fcomp(beta_index,i));
        GP(i) = sum(TP_fcomp(gamma_index,i));
        EP(i) = sum(TP_fcomp(epsilon_index,i));
    end

    norm_ThP = (ThP/mean(ThP(norm_index)))*100;
    norm_AP = (AP/mean(AP(norm_index)))*100;
    norm_BP = (BP/mean(BP(norm_index)))*100;
    norm_GP = (GP/mean(GP(norm_index)))*100;
    norm_EP = (EP/mean(EP(norm_index)))*100;

%% Save
Results.trials = trials;
Results.trialsSpec = trialsSpec;
Results.trial_LFP = trial_LFP;
Results.TP = ThP;
Results.AP = AP;
Results.BP = BP;
Results.GP = GP;
Results.EP = EP;
Results.norm_TP = norm_ThP;
Results.norm_AP = norm_AP;
Results.norm_GP = norm_GP;
Results.norm_BP = norm_BP;
Results.norm_GP = norm_GP;
Results.norm_EP = norm_EP;

Results.RT_hit = RT_hit;
Results.RT2_hit = RT2_hit;
Results.RT_FA = RT_FA;
Results.RT2_FA = RT2_FA;
Results.RT_CR = RT_CR;
Results.RT2_CR = RT2_CR;


Results.time = time;
Results.freq = mor_hz;
save([string(mouseID)+' '+LFP_Channel+' '+graph+' '+t0+' '+spec_type+'power_GNG_multifileAnalysisV2'+'.mat'], 'Results');
    
    % figure;
    %     plot(time-pre,movmean(norm_ThP,smf));
    %     xline(0, '--','LineWidth',1,'Label','NOSE POKE ENTRY', 'LabelHorizontalAlignment', 'Left', 'Color', 'black');
    %     xlabel('Time (seconds)');
    %     title((string(mouseID)) + ' ' + LFP_Channel + ' GNG Normalized Theta Power (6-10 Hz)', graph + ' (' + num2str(size(trial_LFP,2)) + ' trials)');
    % 
    % figure;
    %     plot(time-pre,movmean(norm_AP,smf));
    %     xline(0, '--','LineWidth',1,'Label','NOSE POKE ENTRY', 'LabelHorizontalAlignment', 'Left', 'Color', 'black');
    %     xlabel('Time (seconds)');
    %     title((string(mouseID)) + ' ' + LFP_Channel + ' GNG Normalized Alpha Power (10-15 Hz)', graph + ' (' + num2str(size(trial_LFP,2)) + ' trials)');
    % 
    % figure;
    %     plot(time-pre,movmean(norm_BP,smf));
    %     xline(0, '--','LineWidth',1,'Label','NOSE POKE ENTRY', 'LabelHorizontalAlignment', 'Left', 'Color', 'black');
    %     xlabel('Time (seconds)');
    %     title((string(mouseID)) + ' ' + LFP_Channel + ' GNG Normalized Beta Power (15-30 Hz)', graph + ' (' + num2str(size(trial_LFP,2)) + ' trials)');
    % 
    % figure;
    %     plot(time-pre,movmean(norm_GP,smf));
    %     xline(0, '--','LineWidth',1,'Label','NOSE POKE ENTRY', 'LabelHorizontalAlignment', 'Left', 'Color', 'black');
    %     xlabel('Time (seconds)');
    %     title((string(mouseID)) + ' ' + LFP_Channel + ' GNG Normalized Gamma Power (30-90 Hz)', graph + ' (' + num2str(size(trial_LFP,2)) + ' trials)');
    % 
    % figure;
    %     plot(time-pre,movmean(norm_EP,smf));
    %     xline(0, '--','LineWidth',1,'Label','NOSE POKE ENTRY', 'LabelHorizontalAlignment', 'Left', 'Color', 'black');
    %     xlabel('Time (seconds)');
    %     title((string(mouseID)) + ' ' + LFP_Channel + ' GNG Normalized Epsilon Power (90-150 Hz)', graph + ' (' + num2str(size(trial_LFP,2)) + ' trials)');
%end

clear trialsSpec Results
end

%% Extract/clean up EEG and resp for PAC

if PAC == 'Y'
    OFF_Phase = [];
    OFF_Amp = [];

    hits_LFP = [];
    misses_LFP = [];
    correctRejects_LFP = [];
    falseAlarms_LFP = [];

    hits_resp = [];
    misses_resp = [];
    correctRejects_resp = [];
    falseAlarms_resp = [];

    h = 1;
    m = 1;
    f = 1;
    c = 1;
    
    for i = 1:size(trials,1);
        if trials.Response(i) == "hit" %&& trials.delta(i) >= pre_pac
            %eval('hits_LFP{h} = data(trials.start(i)*fs : trials.end(i)*fs);')
            hits_LFP{h} = data((trials.end(i) - pre_pac)*fs : trials.end(i)*fs);
            hits_resp{h} = resp_data((trials.end(i) - pre_resp)*fs : trials.end(i)*fs);
            hits_delta(h) = trials.delta(i);
            h = h + 1;
        elseif trials.Response(i) == "correctReject" %&& trials.delta(i) >= pre_pac
            correctRejects_LFP{c} = data((trials.end(i) - pre_pac)*fs : trials.end(i)*fs);
            correctRejects_resp{c} = resp_data((trials.end(i) - pre_resp)*fs : trials.end(i)*fs);
            correctRejects_delta(c) = trials.delta(i);
            c = c + 1;
        elseif trials.Response(i) == "miss"%&& trials.delta(i) >= pre_pac
            misses_LFP{m} = data((trials.end(i) - pre_pac)*fs : trials.end(i)*fs);
            misses_resp{m} = resp_data((trials.end(i) - pre_resp)*fs : trials.end(i)*fs);
            misses_delta(m) = trials.delta(i);
            m = m + 1;
        elseif trials.Response(i) == "falseAlarm"%&& trials.delta(i) >= pre_pac
            falseAlarms_LFP{f} = data((trials.end(i) - pre_pac)*fs : trials.end(i)*fs);
            falseAlarms_resp{f} = resp_data((trials.end(i) - pre_resp)*fs : trials.end(i)*fs);
            falseAlarms_delta(f) = trials.delta(i);
            f = f + 1;
        end
    end

   

    trial_LFP_PAC = [];

    if graph == "HITS"
        trial_LFP_PAC = hits_LFP;
        trial_delta_PAC = hits_delta;
    elseif graph == "CORRECT REJECTS"
        trial_LFP_PAC = correctRejects_LFP;
        trial_delta_PAC = correctRejects_delta;
    elseif graph == "MISSES"
        trial_LFP_PAC = misses_LFP;
        trial_delta_PAC = misses_delta;
    elseif graph == "FALSE ALARMS"
        trial_LFP_PAC = falseAlarms_LFP;
        trial_delta_PAC = falseAlarms_delta;
    end

    if remove_artifacts == 'Y'
       trial_LFP_PAC(:,artifacts) = [];
       trial_delta_PAC(artifacts) = [];
    end

%% remove overlapping trials

OL_trial_idx = find(trial_delta_PAC < pre_pac*2);
trial_LFP_PAC(:,OL_trial_idx) = [];


    freq1 = phase_freq(1):dfreq:phase_freq(2);
    freq2 = amp_freq(1):dfreq:amp_freq(2);
    scales_a = frq2scale(freq1,['cmor' num2str(fb) '-' num2str(fc)], 1/fs);
    scales_b = frq2scale(freq2,['cmor' num2str(fb) '-' num2str(fc)], 1/fs);

    filt_resp = [];

    if respiration == "Y"
        if graph == "HITS"
            filt_resp = hits_resp;
        elseif graph == "CORRECT REJECTS"
            filt_resp = correctRejects_resp;
        elseif graph == "MISSES"
            filt_resp = misses_resp;
        elseif graph == "FALSE ALARMS"
            filt_resp = falseAlarms_resp;
        end

         if remove_artifacts == 'Y'
            filt_resp(:,artifacts) = [];
         end
     filt_resp(:,OL_trial_idx) = [];
        for x = 1:length(filt_resp)
            filt_resp{x} = detrend(filt_resp{x});  
            filt_resp{x} = eegfilt(filt_resp{x}', fs, 1 , 25, 0, 500); % from 500 %% FIX THIS! This is the problematic line, the last number needs to be changed!
        end
    end

%% PAC ANALYSIS

    if size(trial_LFP_PAC,2) > 1
        for x = 1:size(trial_LFP_PAC,2)
            %PAC
            % phase analysis
            [p_cwt,mor_hz1] = cwt(trial_LFP_PAC{x},scales_a,['cmor' num2str(fb) '-' num2str(fc)], 1/fs); % Wavelet transform 
            OFF_Phase = cat(2,OFF_Phase,angle(p_cwt));
            clear p_cwt
    
            % amp data
            [a_cwt,mor_hz2] = cwt(trial_LFP_PAC{x},scales_b,['cmor' num2str(fb) '-' num2str(fc)], 1/fs); % Wavelet transform 
            OFF_Amp = cat(2,OFF_Amp, abs(a_cwt).^2);
            clear a_cwt
        end
    end
    
    nbin = 18; % % we are breaking 0-360o in 18 bins, ie, each bin has 20o
    position = zeros(1,nbin); % this variable will get the beginning (not the center) of each bin(in rads)
    winsize = 2*pi/nbin;
    
    for j = 1:nbin
        position(j) = -pi+(j-1)*winsize;
    end
    
    OFF_MI=NaN(size(OFF_Amp,1),size(OFF_Phase,1));  
    
    for x = 1:size(OFF_Phase,1)
        I = [];
        for j = 1:nbin
            I{j}(:) = find(OFF_Phase(x,:) <  position(j) + winsize & OFF_Phase(x,:) >=  position(j));
        end

        for y = 1:size(OFF_Amp,1)
            % now we compute the mean amplitude in each phase:
            MeanAmp = zeros(1,nbin);
            for i = 1:nbin
                MeanAmp(i) = nanmean(OFF_Amp(y,I{i}));
            end
            ind = []; 
            ind = find(isnan(MeanAmp));
            MeanAmp(ind) = 0;
                
            % The center of each bin (for plotting purposes) is
            % position+winsize/2
            % Quantify the amount of amp modulation by means of a
            % normalized entropy index (Tort et al PNAS 2008):
            
            OFF_MI(y,x) = inv_entropy_no_save(MeanAmp');
            %MI(x,y,p)=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);
        end
    end

    figure;
        imagesc(mor_hz1, sort(mor_hz2,'ascend'), ((OFF_MI)));
        eval(['colormap ' cmap ';']);
        set(gca,'YDir','normal');
        %caxis ([0 10e-3]); %caxis([0 30e-3]);
        title((string(mouseID)) + ' ' + LFP_Channel + ' PAC Comodulogram', graph + ' (' + num2str(size(trial_LFP_PAC,2)) + ' trials)');

    if respiration == 'Y'
    
        params.tapers = [2 3]; % # of tapers and bandwith = (2*tapers) -1)
        params.fpass = [0 20]; %frequency band of interest
        params.Fs = fs;
        params.pad = 0;
             
        for x = 1:size(filt_resp,2)
            [resp_psd{x},f] = mtspectrumc(filt_resp{x},params);
            %[resp_psd(x,:),f] = pwelch(analysis_resp(x,:),fs);
        end
        
        for i = 1:length(resp_psd)
            x1 = 1:1:size(resp_psd{i},1);  
            v1 = (resp_psd{i})';
            xq1 = 1:1:250;  %Define num of points for all spectra (more than max... hopefully...
            resp_psd_new(:,i) = (interp1(x1,v1,xq1, 'linear', 'extrap'))';    %Interpolate the function at the query points.
            clear x1 v1 
        end
        
        ave_resp_psd = mean(resp_psd_new,2);
        freq = params.fpass(2)/250:params.fpass(2)/250:params.fpass(2);

        figure;
            plot(freq, ave_resp_psd, ":", 'LineWidth', 2);
            xlim(phase_freq)

%% Save
Results.filt_resp = filt_resp;
Results.resp_psd = resp_psd;

Results.freq = freq;
save([string(mouseID)+' '+LFP_Channel+' '+graph+'_GNG_RESP_multifileAnalysis'+'.mat'], 'Results');

clear Results
    end
    
    Results.MI = OFF_MI;
    Results.phaseHz = mor_hz1;
    Results.ampHz = mor_hz2;
    save([string(mouseID)+' '+LFP_Channel+' '+graph+'_GNG_PAC_multifileAnalysis'+'.mat'], 'Results');
end


