clear

toolboxFolder = ['/Users/emilia/Documents/MATLAB']; %path to the DOTHUB toolbox
addpath(genpath(toolboxFolder)); 

experimentPath = ['/Users/emilia/Documents/LUMO data']; %path where the data is saved
addpath(genpath(experimentPath));

LUMODirName = 'filename.LUMO';

[nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs(LUMODirName);
    
%DOTHUB_dataQualityCheck(nirsFileName);

dod = hmrIntensity2OD(nirs.d);
SD3D = enPruneChannels(nirs.d,nirs.SD3D,ones(size(nirs.t)),[0 1e11],12,[0 100],0); 

SD3D = DOTHUB_balanceMeasListAct(SD3D);

SD2D = nirs.SD; 
SD2D.MeasListAct = SD3D.MeasListAct;

dod = hmrBandpassFilt(dod,nirs.t,0,0.5);
dc = hmrOD2Conc(dod,SD3D,[6 6]);
dc = dc*1e6;

event_code = 'a'
col = find(strcmp(nirs.CondNames, event_code));
events = nirs.s(:,col);
        
[Avg,AvgStd,HRF] = hmrBlockAvg(dc,events,nirs.t,[-5 15]);

dodRecon = DOTHUB_hmrConc2OD(Avg/1e6,SD3D,[6 6]);
tDOD = HRF;

[pathstr, name, ~] = fileparts(nirsFileName);
ds = datestr(now,'yyyymmDDHHMMSS');
preproFileName = fullfile(pathstr,[name '.prepro']);
logData(1,:) = {'Created on: '; ds};
logData(2,:) = {'Derived from data: ', nirsFileName};
logData(3,:) = {'Pre-processed using:', mfilename('fullpath')};
[prepro, preproFileName] = DOTHUB_writePREPRO(preproFileName,logData,dodRecon,tDOD,SD3D,nirs.s,Avg,AvgStd,HRF,nirs.CondNames,SD2D);

condition = 1;
y = squeeze(prepro.dcAvg(:,:,:,condition));
figure(2);
DOTHUB_LUMOplotArray(y,prepro.tHRF,prepro.SD2D, distRange=[20 40]);
