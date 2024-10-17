
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This script runs both time-frequency and phase-amp coupling analyses   %
%  for Intan mice with or without thermistor data.                        %
%  Trials can be organized by reaction time.                              %
%  For PAC and respiratory trace analyses, trials where the hold time is  % 
%  >= 3 seconds are used for the purpose of (hopefully) generating more   %
%  reliable respiratory traces.                                           %
%  If PAC = 'Y', TF plots will be generated only for trials included in   %
%  PAC/resp analyses.                                                     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%% Data Entry

% enter filename  
filelist = 'Therm13multifile';
mouseID = 'Therm13';

%% Behavioral response of interest

graph = "CORRECT"; % PREMATURE, CORRECT, OMISSION

%% LFP channel to analyze
% All channels used %
label = {'CA1', 'BF', 'OB', 'Prl', 'RS', 'AI', 'V1'}'; % CHECK TO ENSURE ORDER IS CORRECT...

%% T=0, Time-locking trials

lever_cond = "LEVER OFF"; % LEVER ON, LEVER OFF, CUE ON

%% Spectral Params & Bin Size

% For T=0 --> Lever OFF

pre = 5; % Seconds prior to lever press/cue on/lever release for time-frequency analyses
post = 0; % Seconds after lever press/cue on/lever release for time-frequency analyses

% Attention: pre = 2, post = 0 (two seconds before lever release)
% Post-attention: pre = -2, post = 4 (2-4 seconds after lever release)

smf = 50; % smoothing factor # of pts

%% Connectivity analysis method



%% Analyze by reaction time?

RT_analysis = "ALL"; % ALL, TOP PERCENTILE, BOTTOM PERCENTILE
RT_PAC = "ALL"; % ALL, TOP PERCENTILE, BOTTOM PERCENTILE
percentile = 25;

%% Filter & remove trial records with artifacts

zScore = 'Y'; % Y or N
notch_filter = 'Y'; % Y or N
Jade = 'Y'; % Y or N
remove_artifacts = 'Y'; % Y or N
show_outliers = 'N'; % Y or N
percent_sat = 5; % percent of total trial length at or above max physiological voltage
%max_v = 350; % max physiological voltage for channel in uV % Now defined below

%% Downsample data?

down_sample = 'Y';
ds = 4; % downsample data (Fs/ds)

%% Plotting parameters

cmap = 'jet'; %color map for spectral plot 

%% Import data

[num, txt] = xlsread([filelist,'.xlsx']);
filesLFP = txt(:,1);
filesTTL = txt(:,2);

reward_channel = 4; % head in reward port 
cue_channel = 5; % cue light on
lever_channel = 3; % Lever depressed
resp_channel = 6; % respiration trace

data = [];
reward_times = [];
cue_times = [];
lever_times = [];
resp_data = [];

endtime = 0;

for y = 1:length(filesLFP)
    fileLFP = load([filesLFP{y}]);
    fileTTL = load([filesTTL{y}]);

    eval(['pdata = fileLFP.LFP;']); % EEG record (mV)
    data = cat(2,data,pdata);
    %eval(['fs = fileLFP.info_lfp.header.sampleRate;']); % EEG sampling frequency Hz (2000)
    fs = 2000; % Just in case there are still sampling rate issues in the OEtoMat script
    eval(['preward_times = fileTTL.' filesTTL{y} '_Ch' num2str(reward_channel) '.times;']); % time (secs)
    reward_times = cat(1,reward_times,(preward_times + endtime));
    eval(['pcue_times = fileTTL.' filesTTL{y} '_Ch' num2str(cue_channel) '.times;']); % binary
    cue_times = cat(1,cue_times,(pcue_times + endtime));
    eval(['plever_times = fileTTL.' filesTTL{y} '_Ch' num2str(lever_channel) '.times;']); % time (secs)
    lever_times = cat(1,lever_times,(plever_times + endtime));
    eval(['presp_data = fileTTL.' filesTTL{y} '_Ch' num2str(resp_channel) '.values*1000;']); % time (secs)
    resp_data = cat(1,resp_data,presp_data);
    endtime=length(data)*(1/fs);
end

clear pdata preward_times pcue_times plever_times presp_data

%% Clean up data

data = detrend(data');

% data = (data');

max_v = 3*rms(data); % Change max voltage value for artifact filtration depending on the channel

% 60/120/180 Hz notch filters
if notch_filter == "Y"
    d1 = designfilt('bandstopiir','FilterOrder',10, ...
       'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5, ...
       'DesignMethod','butter','SampleRate',fs);
    
    d2 = designfilt('bandstopiir','FilterOrder',10, ...
       'HalfPowerFrequency1',119.5,'HalfPowerFrequency2',120.5, ...
       'DesignMethod','butter','SampleRate',fs);
    
    % d3 = designfilt('bandstopiir','FilterOrder',10, ...
    %    'HalfPowerFrequency1',179.5,'HalfPowerFrequency2',180.5, ...
    %    'DesignMethod','butter','SampleRate',fs);
    
    data = filtfilt(d1,data);
    data = filtfilt(d2,data);
    % data = filtfilt(d3,data);

end    

% Downsample data and time

if down_sample == 'Y'
    data = resample(data,fs/ds,fs);
    resp_data = resample(resp_data,fs/ds,fs);
    fs = fs/ds; % convert fs to match
end

data = data';

% if zScore == 'Y'
%     for i = 1:size(data,1)
%         data(i,:) = zscore(data(i,:));
%     end
% end

%% Make table of trial times (Lever Press Times) for all trials 

%Lever Press
Lever_count = 1;
for i = 1:2:length(lever_times)-1
    LeverPress_time(Lever_count) = lever_times(i);
    LeverRelease_time(Lever_count) = lever_times(i + 1);
    Lever_count = Lever_count + 1;
end

% Manual debounce of data
trial_count = 1;
leverON_time = [];
leverOFF_time = [];

for x = 1:length(LeverPress_time)
    if x == 1
        leverON_time(trial_count) = LeverPress_time(x);
        a = find(LeverRelease_time > LeverPress_time(x)+0.05,1);
        leverOFF_time(trial_count)= LeverRelease_time(a);
        trial_count = trial_count + 1;
    elseif LeverPress_time(x) - LeverPress_time(x - 1) > 0.05 && LeverPress_time(x) < endtime - 5 % last part added to prevent issue of lever-pressing when task ends;
        leverON_time(trial_count) = LeverPress_time(x);                                            % can change the last number (seconds) if desired
        a = find(LeverRelease_time > LeverPress_time(x)+0.05,1); %%%                               % has to be changed from 5 to 25 for Therm19 IF and only IF Therm19 is the last mouse in the analysis! :(
        leverOFF_time(trial_count) = LeverRelease_time(a);                                                    
        trial_count = trial_count + 1;
    end
end

% Cue Light
Cue_count = 1;
for i = 1:2:length(cue_times) - 1
    CueOn_time(Cue_count) = cue_times(i);
    CueOff_time(Cue_count) = cue_times(i + 1);
    Cue_count = Cue_count + 1;
end

% Rewards
Reward_count = 1;
for i = 1:2:length(reward_times) - 1
    if reward_times(i + 1) - reward_times(i) > 0.05
        RewardEntry_time(Reward_count) = reward_times(i);
        RewardExit_time(Reward_count) = reward_times(i + 1);
        Reward_count = Reward_count + 1;
    end
end

%% Determine index of each trial type

Premature_trial = [];
Omission_trial = [];
Correct_trial = [];

% switch lever_cond
%     case "LEVER ON"
        for x = 1:length(leverON_time)
            if leverON_time(x) < CueOn_time(end)
                a = find(CueOn_time > leverON_time(x),1);
                %if CueOn_time(a) < leverOFF_time(x) && CueOn_time(a) - leverOFF_time(x) < 5 % Correct
                if CueOn_time(a) < leverOFF_time(x) && leverOFF_time(x) - CueOn_time(a) < 5 % Correct
                    Correct_trial(end+1) = leverON_time(x);
                % elseif CueOn_time(a) < leverOFF_time(x) && CueOn_time(a) - leverOFF_time(x) > 5 % Omission
                elseif CueOn_time(a) < leverOFF_time(x) && leverOFF_time(x) - CueOn_time(a) > 5 % Omission
                    Omission_trial(end+1) = leverON_time(x);
                else 
                    Premature_trial(end+1) =  leverON_time(x);
                end
            end
        end
            
%     case "LEVER OFF"
%         for x = 1:length(leverON_time)
%             if leverON_time(x) < CueOn_time(end)
%                 a = find(CueOn_time > leverON_time(x),1);
%                 if CueOn_time(a) < leverOFF_time(x) && CueOn_time(a) - leverOFF_time(x) < 5 % Correct
%                     Correct_trial(end+1) = leverOFF_time(x); % was leverON_time(x)?
%                 elseif CueOn_time(a) < leverOFF_time(x) && CueOn_time(a)-leverOFF_time(x) > 5 % Omission
%                     Omission_trial(end+1) = leverOFF_time(x);
%                 else 
%                 Premature_trial(end+1) = leverOFF_time(x); % Premature           
%                 end
%             end
%         end
% end

%% Derive RTs

RT = [];
RT_idx = [];
RT_idx_keep = [];
RT_keep = [];
% RT_idx_omiss = [];
% RT_omiss = [];
c = 1;
% d = 1;

for i = 1:length(CueOn_time)-1 % -1 for omissions
    RT_idx(i) = find(leverOFF_time > CueOn_time(i),1);
    RT(i) = leverOFF_time(RT_idx(i)) - CueOn_time(i); % Lever Release

    if RT(i) < 5
        RT_idx_keep(c) = RT_idx(i);
        RT_keep(c) = RT(i);
        c = c + 1;
    end

    % if RT(i) > 5 && RT(i) < 10
    %     RT_idx_omiss(d) = RT_idx(i);
    %     RT_omiss(d) = RT(i);
    %     d = d + 1;
    % end
end

if graph == 'CORRECT'
    RT_idx = RT_idx_keep;
    RT = RT_keep;
else graph == 'OMISSION'
    RT_idx = RT_idx_omiss;
    RT = RT_omiss;
end

for i = 1:length(CueOn_time) 
    RT_idx(i) = find(leverOFF_time > CueOn_time(i),1); 
end

%% Extract EEG 

% For each behavioral type (Correct, Premature, Omission)
switch lever_cond
    case "LEVER ON"
        switch graph 
            case "CORRECT"
                trial_time = leverON_time(RT_idx);
            
            case "PREMATURE"
                trial_time = Premature_trial;
            
            case "OMISSION"
                trial_time = Omission_trial;
        end

    case "LEVER OFF"
        switch graph 
            case "CORRECT"
                trial_time = leverOFF_time(RT_idx);
            
            case "PREMATURE"
                trial_time = Premature_trial;
            
            case "OMISSION"
                trial_time = Omission_trial;
        end
        
    case "CUE ON"
        switch graph 
            case "CORRECT"
                trial_time = CueOn_time;
            
            case "PREMATURE"
                trial_time = Premature_trial;
            
            case "OMISSION"
                trial_time = Omission_trial;
        end
                
end

% Dividing EEG data into trials

trial_LFP = [];

for i = 1:length(trial_time)
    trial_LFP(:,:,i) = data(:,(trial_time(i) - pre)*fs :(trial_time(i) + post)*fs); % Data time-locked to lever press
end

clear data

%% Signal and artifact filtering

if remove_artifacts == 'Y'

     % Channel indexing (for removal of "bad" channels during artifact filtering while retaining structure of trial_LFP)
     % 1-CA1, 2-BF, 3-OB, 4-Prl, 5-RS, 6-AI, 7-V1

    if mouseID == 'Therm13' | mouseID == 'Therm20'
        ch_idx = [2 3 4 5 6 7]; % CA1 removed
    elseif mouseID == 'Therm16'  % Plots with RS come out flat, but granger matrix shows RS activity?
        ch_idx = [1 2 3 4 6 7]; % RS removed
    else mouseID == 'Therm19'
        ch_idx = [1 2 3 4 5 6 7]; % none removed
    end

    for x = 1:length(ch_idx)
        y = ch_idx(x);
    end

    % JADE signal processing

    if Jade == 'Y'
        for i = 1:size(trial_LFP,3)
            [A,trial_LFP(y,:,i)] = jade(trial_LFP(y,:,i));
        end
    end

    % Identifying and removing trials with artifacts in any channel from future analyses

    for i = 1:size(trial_LFP,1)

        a = 1;
        sf = []; 
        ent = []; 
        max_pts = []; 
        maxx = []; 
        activity = []; 
        mobility = []; 
        complexity = [];

    
        for j = 1:size(trial_LFP,3)
            sf(j) = sum(spectralFlatness(trial_LFP(y,:,j)',fs));
            [activity(j), mobility(j), complexity(j)] = hjorth(trial_LFP(y,:,j));
            ent(j) = mean(wentropy(trial_LFP(y,:,j)));
            %FNN(j) = f_fnn(analysis_LFP(y,;,j),2,1,15,2);
            b = find(abs(trial_LFP(y,:,j))>=max_v(y)); % Not relevant here bc max_v was calculated prior to Jade (?)
            max_pts(j) = length(b)/size(trial_LFP,3) * 100;
            if length(b)/size(trial_LFP,3)>percent_sat/100 % percentage of saturation above which trial will be omitted from spectral analysis
                maxx(a) = j;
                a = a + 1;
            end

        end
    
        if show_outliers == 'Y'
            figure;
            plot(sf, 'b.');
            title('Spectral Flatness');
            
            figure;
            plot(activity, 'b.');
            title('Hjorth Activity');
            
            figure;
            plot(mobility, 'b.');
            title('Hjorth Mobility');
            
            figure;
            plot(complexity, 'b.');
            title('Hjorth Complexity');
            
            figure;
            plot(ent, 'b.');
            title('Shannon Entropy');
    
            figure;
            plot(max_pts, 'b.');
            title('Percent of Data Points Above Max Voltage');
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
            disp('Spectral Flatness Outliers:');
            disp(err_sf);
            disp('Activity Outliers:');
            disp(err_activity);
            disp('Complexity Outliers:');
            disp(err_complexity);
            disp('Entropy Outliers:');
            disp(err_ent);
            disp('Max Voltage Outliers:');
            disp(maxx);
        end

        %bad_trials = [err_sf, err_activity, err_mobility, err_complexity, err_ent, maxx]; 
        bad_trials = [err_sf, err_activity, err_mobility, err_complexity, err_ent]; % removed maxx
        artifacts = unique(bad_trials);
        trial_LFP(:,:,artifacts) = []; % removes any trial that is flagged from all channels, I think??

    end

    % Z-score

    if zScore == 'Y'
        for i = 1:size(trial_LFP,3)
            trial_LFP(y,:,i) = zscore(trial_LFP(y,:,i));
        end

    end
end

%%  Analyze by RT

RT = [];

for i = 1:length(CueOn_time)
    RT(i) = leverOFF_time(RT_idx(i)) - CueOn_time(i); % Lever Release
end

RT(:,artifacts) = []; % Removes trials with artifacts from RT analyses

switch RT_analysis
    case "ALL"
        analysis_LFP = trial_LFP;
    case "TOP PERCENTILE"
        [a, idx] = mink(RT, round(length(RT)*(percentile/100))); % lever release
        % a = find(CueOn_time > leverON_time(:,idx(x)), "first");
        for x = 1:length(idx)
            analysis_LFP(:,x) = trial_LFP(:,idx(x));
        end
    case "BOTTOM PERCENTILE"
        [a, idx] = maxk(RT, round(length(RT)*(percentile/100))); % lever release
        for x = 1:length(idx)
            analysis_LFP(:,x) = trial_LFP(:,idx(x));
        end
end  
    
%% Connectivity Analysis 
clear pdata presp_data
data = []; dtime = 1/fs:1/fs:size(trial_LFP,2)/fs;

%% Structure data for field trip analysis of connectivity & Zscore data

data = []; 
dtime = 1/fs:1/fs:size(trial_LFP,2)/fs;
data.label = label; 
data.fsample = fs;

for x = 1:size(trial_LFP,3)
    for y = 1:size(trial_LFP,1)
        data.trial{x}(y,:) = (trial_LFP(y,:,x));
        data.time{x}(1,:) = dtime;
    end
end

switch cfg.method
    case  'mi'
        mi = ft_connectivityanalysis(cfg, data);
        figure; 
        hold on;
        % imagesc(abs(mi.mi)^2); %% NOT SURE IF THIS IS PROPER WAY TO TRANSFORM COMPLEX VALUES
        imagesc(real(mi.mi)); %%
        colorbar('fontsize', 20);
        set(gca,'fontsize', 20, 'Xtick', 10:10:40, 'Ytick', 10:10:40, 'tickdir', 'out');
        xticks(1:numel(label));
        xticklabels(label);
        yticks(1:numel(label));
        yticklabels(label);
     
    case 'granger'
        
        cfg = [];
        cfg.order = 2;
        mdata = ft_mvaranalysis(cfg, data);
        mdata.offset = 0;
        cfg = [];
        cfg.method = 'mvar';
        cfg.foi = [0:200];
        mfreq = ft_freqanalysis(cfg, mdata);
        %calculate the fourier coefficients (non-parametric derivation of power)
        cfg = []; 
        % cfg.method = 'wavelet';
        % cfg.foilim = [1 150]; %frequency band of interest
        % cfg.toi    = 0.5; % vector 1 x numtoi, the times on which the analysis windows should be centered (in seconds)
        % cfg.width  = 13; % 'width', or number of cycles, of the wavelet (default = 7)
        % cfg.gwidth = 3; % determines the length of the used wavelets in standard deviations of the implicit Gaussian kernel and should be chosen >= 3; (default = 3)
        cfg.method = 'mtmfft';
        cfg.taper = 'dpss';
        cfg.output = 'fourier';
        cfg.tapsmofrq = 4;
        cfg.foilim = [0 100];
        cfg.pad = 'nextpow2';
        freq = ft_freqanalysis(cfg, data);
        cfg = [];
        cfg.method = 'coh';
        cfg.complex = 'abs';
        coh1 = ft_connectivityanalysis(cfg, freq);
        csd = ft_checkdata(freq, 'cmbrepresentation', 'fullfast');
        cfg = [];
        cfg.order = 2;
        mdata = ft_mvaranalysis(cfg, data);
        mdata.offset = 0;
        cfg = [];
        cfg.method = 'mvar';
        cfg.foi = [0:120];
        mfreq = ft_freqanalysis(cfg, mdata);
        cfg = [];
        cfg.method = 'coh';
        cohp = ft_connectivityanalysis(cfg, mfreq);
        cfg = [];
        cfg.method = 'granger';
        cfg.granger.sfmethod = 'bivariate';
        g1 = ft_connectivityanalysis(cfg, csd);
        cfg = [];
        cfg.method = 'granger';
        gp = ft_connectivityanalysis(cfg, mfreq);
        
       
        figure;plot(g1.freq,squeeze(g1.grangerspctrm(33,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(34,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('Prl->AI','AI->Prl');
        figure;plot(g1.freq,squeeze(g1.grangerspctrm(15,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(16,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('BF->PrL','PrL->BF');
        figure;plot(g1.freq,squeeze(g1.grangerspctrm(19,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(20,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('BF->AI','AI->BF');

        %if mouseID ~= 'Therm13' | mouseID ~= 'Therm20'
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(1,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(2,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('CA1->BF','BF->CA1');
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(3,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(4,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('CA1->OB','OB->CA1');
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(5,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(6,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('CA1->PrL','PrL->CA1');
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(7,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(8,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('CA1->V1','V1->CA1');
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(9,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(10,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('CA1->AI','AI->CA1');
           % if mouseID ~= 'Therm16'
                figure;plot(g1.freq,squeeze(g1.grangerspctrm(11,:)));hold on
                plot(g1.freq,squeeze(g1.grangerspctrm(12,:)),'r');ylim([0 0.25]);xlim([0 100]);
                title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('CA1->RS','RS->CA1');
           % end
        %end

        figure;plot(g1.freq,squeeze(g1.grangerspctrm(13,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(14,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('BF->OB','OB->BF');
        figure;plot(g1.freq,squeeze(g1.grangerspctrm(17,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(18,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('BF->V1','V1->BF');

        %if mouseID ~= 'Therm16'
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(21,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(22,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('RS->BF','BF->RS');
        %end

        figure;plot(g1.freq,squeeze(g1.grangerspctrm(23,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(24,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('OB->PrL','PrL->OB');
        figure;plot(g1.freq,squeeze(g1.grangerspctrm(25,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(26,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('OB->V1','V1->OB');
        figure;plot(g1.freq,squeeze(g1.grangerspctrm(27,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(28,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('OB->AI','AI->OB');

        %if mouseID ~= 'Therm16'
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(29,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(30,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('OB->RS','RS->OB');
        %end

        figure;plot(g1.freq,squeeze(g1.grangerspctrm(31,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(32,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('PrL->V1','V1->PrL');

        %if mouseID ~= 'Therm16'
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(35,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(36,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('PrL->RS','RS->PrL');
        %end

        figure;plot(g1.freq,squeeze(g1.grangerspctrm(37,:)));hold on
        plot(g1.freq,squeeze(g1.grangerspctrm(38,:)),'r');ylim([0 0.25]);xlim([0 100]);
        title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('V1->AI','AI->V1');

        % if mouseID ~= 'Therm16'
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(39,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(40,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('V1->RS','RS->V1');
            figure;plot(g1.freq,squeeze(g1.grangerspctrm(41,:)));hold on
            plot(g1.freq,squeeze(g1.grangerspctrm(42,:)),'r');ylim([0 0.25]);xlim([0 100]);
            title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('AI->RS','RS->AI');
        % end

        
%% generate plot of mean causality at beta (15-30 Hz) and gamma bands (30-80 Hz) between each channel

granger_matrix = zeros(length(label));
granger_matrix2 = zeros(length(label));
%freq_idx = find(g1.freq > 30 & g1.freq < 50); % original, 40 +/- 10 Hz
freq_idx = find(g1.freq >= 15 & g1.freq <= 30); % Beta frequency index
freq_idx2 = find(g1.freq >= 30 & g1.freq <= 80); % Gamma frequency index
ct = 1;
ct2 = 1;

for x = 1:size(label)
    for y = 1:size(label)
        if x == y
            granger_matrix(x,y) = 0;
            granger_matrix2(x,y) = 0;
        else
            if granger_matrix(y,x) == 0
                granger_matrix(y,x) = mean(g1.grangerspctrm(ct,freq_idx));
                granger_matrix(x,y) = mean(g1.grangerspctrm(ct+1,freq_idx));
                ct = ct + 2;
            end
            
            if granger_matrix2(y,x) == 0
                granger_matrix2(y,x) = mean(g1.grangerspctrm(ct2,freq_idx2));
                granger_matrix2(x,y) = mean(g1.grangerspctrm(ct2+1,freq_idx2));
                ct2 = ct2 + 2;
            end
        end
    end
end


% figure;hold on; % Original
% imagesc(granger_matrix);
% eval(['colormap ' cmap ';']);
% plot([0.5 7.5],[0.5 7.5],'w','linewidth',2);
% plot([0.5 7.5],[0.5 7.5],'w--','linewidth',2);
% clim([0 0.25]);
% axis square; axis tight
% set(gca,'fontsize', 25, 'Xtick', 10:10:40, 'Ytick', 10:10:40, 'tickdir', 'out');
% title('Mean 40Hz Granger','fontsize',25);
% colorbar('fontsize', 25);

figure;hold on;
imagesc(granger_matrix);
eval(['colormap ' cmap ';']);
plot([0.5 7.5],[0.5 7.5],'w','linewidth',2);
plot([0.5 7.5],[0.5 7.5],'w--','linewidth',2);
clim([0 0.25]);
axis square; axis tight
set(gca,'fontsize', 20, 'Xtick', 10:10:40, 'Ytick', 10:10:40, 'tickdir', 'out');
title('Mean 15-30Hz Granger','fontsize',20);
colorbar('fontsize', 20);
xticks(1:numel(label));
xticklabels(label);
yticks(1:numel(label));
yticklabels(label);

figure;hold on;
imagesc(granger_matrix2);
eval(['colormap ' cmap ';']);
plot([0.5 7.5],[0.5 7.5],'w','linewidth',2);
plot([0.5 7.5],[0.5 7.5],'w--','linewidth',2);
clim([0 0.25]);
axis square; axis tight
set(gca,'fontsize', 20, 'Xtick', 10:10:40, 'Ytick', 10:10:40, 'tickdir', 'out');
title('Mean 30-80Hz Granger','fontsize',20); 
colorbar('fontsize', 20);
xticks(1:numel(label));
xticklabels(label);
yticks(1:numel(label));
yticklabels(label);

end
