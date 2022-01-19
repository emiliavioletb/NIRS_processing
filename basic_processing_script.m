%% NIRS PRE-PROCESSING SCRIPT.

% This script pre-processes group-level LUMO data and block averages 
% across conditions. Includes short channel regression and 
% visualisation of event-related responses. 
% Created by Emilia with help from Raul Hernandez and Georgina Leadley.

%% Specify paths of toolboxes & raw data 
clear

toolboxFolder = ['/Users/emilia/Documents/MATLAB']; %path to the DOTHUB toolbox
addpath(genpath(toolboxFolder)); 

experimentPath = ['/Users/emilia/Documents/raw_data']; %path where the raw data is saved
addpath(genpath(experimentPath));

%% Loop through participants, pre-process data and block average  

participant_list = dir('*.LUMO');
plotFlag = 0; % suppress output of quality checks?

for i = 1:length(participant_list)
    
    LUMODirName = [sprintf('%d',i), '.LUMO'];
    
    % Convert the .LUMO files to .nirs
    [nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs(LUMODirName);
    
    % Run data quality checks. This step produces various figures to check
    % the data quality including computing power-spectral density, motion
    % 
    if plotFlag
        DOTHUB_dataQualityCheck(nirsFileName);
    end

    % Converting raw data to optical density. 
    dod = hmrIntensity2OD(nirs.d);

    % Prune channels with SNR threshold, source-detector separation range
    % or data range. 
    SD3D = enPruneChannels(nirs.d,nirs.SD3D,ones(size(nirs.t)),[0 1e11],12,[0 100],0); 

    %Force MeasListAct to be the same across wavelengths
    nWavs            = length(SD3D.Lambda);
    tmp              = reshape(SD3D.MeasListAct,length(SD3D.MeasListAct)/nWavs,nWavs);
    tmp2             = ~any(tmp'==0)';
    SD3D.MeasListAct = repmat(tmp2,nWavs,1);

    %Set SD2D
    SD2D = nirs.SD; 
    SD2D.MeasListAct = SD3D.MeasListAct;

    %Bandpass filter.
    dod = hmrBandpassFilt(dod,nirs.t,0,0.5);

    % Convert to concentration.
    dc = hmrOD2Conc(dod,SD3D,[6 6]);
    dc = dc*1e6;

    %Regress short channels
    dc = DOTHUB_hmrSSRegressionByChannel(dc,SD3D,12,1);
    
    % Find events for the conditions. This step can be repeated if you 
    % have multiple tasks or conditions, all you need to know are the 
    % event codes. 
    event_column = find(strcmp(nirs.CondNames, 'event_a'));
    events       = nirs.s(:,event_column);

    % Deconvolution (an alternative to regression by channel)
    % [dcAvg_ded,dcAvgStd_ded,~] = hmrDeconvHRF_DriftSS(dc,events,nirs.t,SD3D,[],ones(size(nirs.t)),[a b],1,1,[1 1],12,1,0,0);
        
    % Block average
    a = -5; % Start of block average
    b = 65; % End of block average
    [Avg,AvgStd,HRF] = hmrBlockAvg(dc,events,nirs.t,[a b]);
 
    
    % Output data for each loop/participant
    Avg_group(i,:,:,:)    = Avg;
    AvgStd_group(i,:,:,:) = AvgStd;

    % Save each participant block
    filename = [sprintf('%d',i)]
    save([filename, 'blockAvg.mat'], 'Avg');
     save([filename, 'blockAvgStd.mat'],'AvgStd');

end

%% File saving
% Save the HRF
save('HRF.mat', 'HRF');

% Save group-level averaged responses
save('Avg_group.mat', 'Avg_group');
save('AvgStd_group.mat', 'shakeAvgStd_group');

%% Creating figures

% Get source-detector distances 
dist = DOTHUB_getSDdists(SD3D);
filt30mm = dist > 27.5 & dist < 32.5;
filt20mm = dist > 17.5 & dist < 22.5;
filt10mm = dist > 7.5 & dist < 12.5;

t = 60; % Block duration

% 30mm separation 
block_30mm = Avg_group(:,:,:,filt30mm);
Av_channels_VFTsem30mm = squeeze(mean(block_30mm,4));
Av_participants_30mm = squeeze(mean(Av_channels_30mm,1));
Std30 = AvgStd_group(:,:,:,filt30mm);
Av_channels_30mmStd = squeeze(mean(Std30,4));
Av_participants_30mmStd = squeeze(mean(Av_channels_30mmStd,1));
SEM_30mm = Av_participants_30mmStd/sqrt(22);
figure(1);
DOTHUB_plotHRF(HRF,Av_participants_30mm,SEM_30mm,t)

% 20mm separation 
block_20mm = Avg_group(:,:,:,filt20mm);
Av_channels_20mm = squeeze(mean(block_20mm,4));
Av_participants_20mm = squeeze(mean(Av_channels_20mm,1));
Std20 = AvgStd_group(:,:,:,filt20mm);
Av_channels_20mmStd = squeeze(mean(Std20,4));
Av_participants_20mmStd = squeeze(mean(Av_channels_20mmStd,1));
SEM_20mm = Av_participants_20mmStd/sqrt(22);
figure(1);
DOTHUB_plotHRF(HRF,Av_participants_20mm,SEM_20mm,t)

% 10mm separation
block_10mm = Avg_group(:,:,:,filt10mm);
Av_channels_10mm = squeeze(mean(VFTsem10,4));
Av_participants_VFTsem_10mm = squeeze(mean(Av_channels_10mm,1));
Std10 = AvgStd_group(:,:,:,filt10mm);
Av_channels_10mmStd = squeeze(mean(Std10,4));
Av_participants_10mmStd = squeeze(mean(Av_channels_10mmStd,1));
SEM_10mm = Av_participants_10mmStd/sqrt(22);
figure(1);
DOTHUB_plotHRF(HRF,Av_participants_10mm,SEM_10mm,60)
