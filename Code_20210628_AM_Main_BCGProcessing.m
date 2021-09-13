% Signal Processing Pipeline 
% Gating Method: ECG-R peak
% BCG Preprocessing: Bandpass filter + Exponential Moving Average
% BCG Feature Extraction: PDF method

% List of interventions:------------------------------
% BL (Baseline)
% MA (Mental Arithmetic)
%NB (N-back)
% SC (Stroop Color)
% CP (Cold Pressor)

%---------------------------------------------------------
clear all; clc; close all;

% Add Path to library--------------------------------
currentfolder = pwd;
addpath(currentfolder,'Library')
addpath(currentfolder,'Library External')
addpath(currentfolder,'Library GT')

% Load Data-------------------------------------------
[~,~,subjectInterventionXLS] = xlsread('subjectIntervention'); % Excel file containing interventions timings and info on whether or not include each intervention
% load('D:\Project - ONR Stress\Data - UMD\Full Protocol\Full Protocol\Jesse_full_protocol_04272021'); SID = 10;
% load('D:\Project - ONR Stress\Data - UMD\Full Protocol\Full Protocol\Weidi_full_protocol_05202021'); SID = 11; fileIndex = 1; % Change SID and fileIndex for each new subject. fileIndex should be set equl to the line # in the excel file containing current subject info subtracted by 1
load('D:\Project - ONR Stress\Data - UMD\Pilot Data\20210623SK'); SID = 1; fileIndex = 2;
% load('D:\Project - ONR Stress\Data - UMD\Pilot Data\20210625NE'); SID = 2; fileIndex = 3;
% load('D:\Project - ONR Stress\Data - UMD\Pilot Data\20210630MP'); SID = 3; fileIndex = 4;
% load('D:\Project - ONR Stress\Data - UMD\Pilot Data\20210703MG'); SID = 4; fileIndex = 5;
% load('D:\Project - ONR Stress\Data - UMD\Pilot Data\20210714GA'); SID = 5; fileIndex = 6;

% Define Data-----------------------------------------
ECGtemp = data(:,6);    % Biopac ECG
BCG_Wtemp = data(:,2);   % Wrist BCG: head-to-foot
SCG_HFtemp = data(:,5); % Chest SCG: head-to-foot
SCGtemp = data(:,7);    % Chest SCG
BPtemp = data(:,8);  % BP nexfin
PPGtemp = data(:,3); % finger PPG

% Incorporate 250 [ms] time delay in BP nexfin-------
%%% Since fs = 1000 therefore 250 [ms] is equal to 250 samples 
nd = length(ECGtemp);
range = 251:nd-250;
ECGtemp = ECGtemp(range);
BCG_Wtemp = BCG_Wtemp(range);
SCGtemp = SCGtemp(range);
SCG_HFtemp = SCG_HFtemp(range);
BPtemp = BPtemp(range+250);
PPGtemp = PPGtemp(range);

% General Parameters------------------------------
fs = 1000;  % Sampling Frequency in [Hz]   

% Initialization using the Excel File---------------
validEventNum = 0;eventIdx=1;
while ~isnan(subjectInterventionXLS{fileIndex+1,2+(eventIdx-1)*4})
    % An intervention is included for analysis if its *note* in the excel
    % file is set to 'Y'
    if subjectInterventionXLS{fileIndex+1,2+(eventIdx-1)*4+3} == 'Y' 
        validEventNum = validEventNum+1;
    end
    eventIdx = eventIdx+1;
end

% Looping over all listed interventions in the Excel file------------
validEventIdx=1;
eventIdx=1;
%%% Define variables
ECG = [];
BCG_W = [];
SCG = [];
SCG_HF = [];
BP = [];
PPG = [];
timesig = [];
startIdxEntr = [];
startIdxCut = [];
while ~isnan(subjectInterventionXLS{fileIndex+1,2+(eventIdx-1)*4}) %eventIdx<=4 %
    name = ['event' num2str(eventIdx)];
    eventTag = ['Event' num2str(validEventIdx)];
    eventData.(name).eventName = subjectInterventionXLS{fileIndex+1,2+(eventIdx-1)*4};
    eventData.(name).startTime = round(subjectInterventionXLS{fileIndex+1,2+(eventIdx-1)*4+1}*60);
    eventData.(name).endTime = round(subjectInterventionXLS{fileIndex+1,2+(eventIdx-1)*4+2}*60);
    eventData.(name).note = subjectInterventionXLS{fileIndex+1,2+(eventIdx-1)*4+3};
    interventionName = eventData.(name).eventName;
    %%% Each valid intervention will be included in the analysis
    if eventData.(name).note == 'Y'
        startTime = eventData.(name).startTime;
        endTime = eventData.(name).endTime;
        startIdxEntr = [startIdxEntr; startTime*fs]; % Intervention start sample # wrt the entire recorded signal
        startIdxCut = [startIdxCut; length(timesig)+1]; % Intervention start sample # wrt the valid portion of the signal
        intTime = startTime:1/fs:endTime;
        if ~isempty(timesig)
            timesig = [timesig, intTime-startTime+timesig(end)+1/fs]; % stores time signal based on valid interventions
        else
            timesig = intTime-startTime+1/fs;
        end
        dataRange = startTime*fs:endTime*fs; % defines the range for the current intervention
        ECG = [ECG;ECGtemp(dataRange)]; % Stores ECG in valid interventions
        BCG_W = [BCG_W;BCG_Wtemp(dataRange)]; % Stores BCG_W in valid interventions
        SCG = [SCG;SCGtemp(dataRange)]; % Stores SCG in valid interventions
        SCG_HF = [SCG_HF;SCG_HFtemp(dataRange)]; % Stores SCG_HF in valid interventions
        BP = [BP;BPtemp(dataRange)]; % Stores BP in valid interventions
        PPG = [PPG;PPGtemp(dataRange)]; % Stores PPG in valid interventions
        validEventIdx = validEventIdx+1;
    end
    eventIdx = eventIdx+1;
end

% Define parameter values------------------------
param = SetParameters(fs,SID);

% ECG Processing and R-Peak detection-----
[u,v]=butter(2,[10 30]/(fs/2));
ECGFiltered = filtfilt(u,v,ECG);
% [~,Rpeak] = findpeaks(ECGFiltered,'minpeakheight',param.ECG.RMinPeakHeight,'minpeakdistance',500,'minpeakwidth',10);
clear temp
[~,temp]=pan_tompkin(ECG,fs,0);    % Extracting ECG R-peaks using Pan-Tompkin method
Rpeak = transpose(temp);
RpkKeep = Rpeak;
RpkKeepInt = Rpeak/fs;

% Plot ECG and detected R-peaks--------------
figure
ax(1) = subplot(6,1,1);
plot(timesig,ECG,'g')
hold on
plot(timesig,ECGFiltered,'k')
plot(RpkKeepInt,ECG(RpkKeep),'ok')
ylabel('ECG')
xticks(RpkKeepInt)
grid on
title(['SID: ',num2str(SID)])
currAxis = gca;
textloc = 1.5; % location of the text on the plot
PlotInterVLine(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,1) % split the plot with vertical line marking the start of each valid intervention

% BP Processingand feature extraction---------
[u,v] = butter(2,20/(fs/2));
BP = filtfilt(u,v,BP);
ft = []; % Variable storing the features
ft = FindBPFeatures2(ft, BP,'BP',param.BP,Rpeak,fs,1,length(Rpeak));

% BP plot-----------------------------------------------
ax(2) = subplot(6,1,2);
plot(timesig,BP,'g')
hold on
plot(ft.DP(:,2)/fs,ft.DP(:,1),'*b')
plot(ft.SP(:,2)/fs,ft.SP(:,1),'*r')
plot(ft.tanDP(:,2)/fs,ft.tanDP(:,1),'ob')
plot(ft.tanSP(:,2)/fs,ft.tanSP(:,1),'or')
ylabel('BP')
xticks(RpkKeepInt)
grid on
currAxis = gca;
PlotInterVLine(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,0) % split the plot with vertical line marking the start of each valid intervention

% PPG Preprocessing and feature extraction--------------
[u,v] = butter(2,20/(fs/2));
PPG = filtfilt(u,v,PPG);
ft = FindBPFeatures2(ft, PPG,'PPG',param.PPG,Rpeak,fs,1,length(Rpeak));

% BCG Preprocessing------------------------------
startIdx = 1;
[u,v]=butter(2,[0.5,10]/(fs/2));
BCGFilt_W = filtfilt(u,v,BCG_W);

% PPG plot--------------------------------------------
ax(3) = subplot(6,1,3);
plot(timesig,PPG,'g')
hold on
plot(ft.PPGFt(:,2)/fs,ft.PPGFt(:,1),'*b')
plot(ft.tanPPGFt(:,2)/fs,ft.tanPPGFt(:,1),'ob')
plot(ft.PPGPk(:,2)/fs,ft.PPGPk(:,1),'*r')
plot(ft.tanPPGPk(:,2)/fs,ft.tanPPGPk(:,1),'or')
ylabel('PPG')
xticks(RpkKeepInt)
grid on
currAxis = gca;
PlotInterVLine(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,0) % split the plot with vertical line marking the start of each valid intervention

% SCG Preprocessing------------------------------
[u,v]=butter(2,[0.5,10]/(fs/2));
% SCG = filtfilt(u,v,SCG);
SCGBut_HF = filtfilt(u,v,SCG_HF);
% GT Filter for SCG and BCG
Hd = createBPF(0.5,40,fs,'fpass1',1,'fpass2', 20, 'kaiser', 'order', 20);
BCGKsr_W = filtfilt(Hd.Numerator,1,BCG_W);
SCGFilt = filtfilt(Hd.Numerator, 1, SCG);
% SCGFilt_HF = filtfilt(Hd.Numerator, 1, SCG_HF);
SCGFilt_HF = SCGBut_HF;

% Plot fft of the SCG---------------------------------
% figure
% calcFFT(SCG,fs);
% figure
% calcFFT(SCGFilt,fs);

% Exclude beats with large motion---------------
MotMult = param.MotMult; % motion artifact exclusion threshold for BCG
MotMultSCG = param.MotMultSCG; % motion artifact exclusion threshold for SCG
[flagBeats, MotLim, flagInd, Rpeak] =  AssessBeats(BCGFilt_W,Rpeak,MotMult); 
[flagBeatsSCG, MotLimSCG, flagIndSCG, RpeakSCG] =  AssessBeats(SCGFilt,RpkKeep,MotMultSCG); 

% Exponential Moving Average on BCG, SCG, and SCG_HF-------------------
BCGMA_W = expMA(BCGFilt_W,BCGFilt_W,Rpeak,startIdx,timesig,fs);
SCGMA = expMA(SCGFilt,SCGFilt,Rpeak,startIdx,timesig,fs);
SCGMA_HF = expMA(SCGFilt_HF,SCGFilt_HF,Rpeak,startIdx,timesig,fs);

% Plot BCG, SCG, and SCG_HF----------------
figure(1)
ax(4) = subplot(6,1,4);
plot(timesig,BCGFilt_W,'g')
hold on
plot(flagInd/fs,BCGFilt_W(flagInd),'r.')
yline(MotLim)
plot(timesig,BCGMA_W,'k')
currAxis = gca;
PlotInterVLine(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,0) % split the plot with vertical line marking the start of each valid intervention
ylabel('BCG')
ax(5) = subplot(6,1,5);
% plot(SCG,'g')
hold on
plot(timesig,SCGFilt,'g')
plot(timesig,SCGMA,'k')
plot(flagIndSCG/fs,SCGFilt(flagIndSCG),'r.')
yline(MotLimSCG)
currAxis = gca;
PlotInterVLine(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,0) % split the plot with vertical line marking the start of each valid intervention
grid on
xticks(RpkKeepInt)
ylabel('SCG')
ax(6) = subplot(6,1,6);
plot(timesig,SCGFilt_HF,'g')
hold on
plot(timesig,SCGMA_HF,'k')
grid on
% plot(timesig,BCGFilt_W,'m')
currAxis = gca;
PlotInterVLine(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,0) % split the plot with vertical line marking the start of each valid intervention
linkaxes(ax,'x')
xlabel('Time [s]')
ylabel('SCG_H_F')
xticks(RpkKeepInt)

% Wrist BCG J-wave extraction-------------------
%%% Phase 1: find j-wave candidates, each beat can have 10 candidates max
ft = FindFeatures(ft, BCGMA_W, 'BCG',param.BCG_A,Rpeak,fs,1,length(Rpeak));
% figure(1)
% subplot(6,1,4)
% plot(ft.J_T/fs,ft.J_A,'o')
J_PEP = ft.J_T-ones(10,1)*transpose(Rpeak(1:end-1));
[~,ind] = max(ft.J_A);
for ii = 1:length(ft.J_A)
J_Phase1(ii) = J_PEP(ind(ii),ii);
end
[f,xi] = ksdensity(J_Phase1);
figure
plot(xi,f)
title(['SID: ',num2str(SID)])
%%% set two ends of PDF equal to zero
[~,indMode] = max(f);
f(abs(xi-xi(indMode))>80) = 0;
hold on
plot(xi,f)
xlabel('Time [ms]')
%%% Phase 2: find the J-wave with highest PDF
[J_Phase2, flag] = FindJ_PDF(J_PEP,ft.J_A,xi,f);

% Extract H, I, K, L waves in wrist BCG------------
ft = FindRestFeat(ft, J_Phase2.Val(:,1:2), BCGMA_W, 'BCG_A',param,Rpeak,fs,1,length(Rpeak));

% Plot extracted BCG Features-------------------
figure(1)
subplot(6,1,4)
hold on
plot((J_Phase2.Val(flag>0,1)+Rpeak(flag>0))/fs,J_Phase2.Val(flag>0,2),'*k')
% plot((J_Phase2.Val(flag==0,1)+Rpeak(flag==0))/fs,J_Phase2.Val(flag==0,2),'*m')
xticks(RpkKeepInt)
% plot(timesig,SCGFilt_HF,'b')
% plot(timesig,SCGBut_HF,'m')
% plot(timesig,BCGKsr_W,'b')
plot(ft.H(:,2)/fs,ft.H(:,1),'om')
plot(ft.I(:,2)/fs,ft.I(:,1),'og')
plot(ft.K(:,2)/fs,ft.K(:,1),'or')
plot(ft.L(:,2)/fs,ft.L(:,1),'ob')
grid on

% Plot J-wave amplitudes on a separate figure-------------
figure
plot((J_Phase2.Val(flag>0,1)+Rpeak(flag>0))/fs,J_Phase2.Val(flag>0,2),'*k')
currAxis = gca;
textloc = 0.8;
PlotInterVLine(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,1) % split the plot with vertical line marking the start of each valid intervention
ylabel('J wave')
title(['SID: ',num2str(SID)])

% Chest SCG Head-to-Foot (Chest BCG) J-wave Extraction--------
%%% Phase 1: find j-wave candidates, each beat can have 10 candidates max
ftSCG_HF = [];
[ftSCG_HF,plotBeat] = FindFeatures(ftSCG_HF, SCGMA_HF, 'BCG',param.SCG,Rpeak,fs,1,length(Rpeak));
figure(1)
subplot(6,1,6)
plot(ftSCG_HF.J_T/fs,ftSCG_HF.J_A,'o')
J_PEP_SCG = ftSCG_HF.J_T-ones(10,1)*transpose(Rpeak(1:end-1));
[~,ind] = max(ftSCG_HF.J_A);
for ii = 1:length(ftSCG_HF.J_A)
J_Phase1_SCG(ii) = J_PEP_SCG(ind(ii),ii);
end
[f,xi] = ksdensity(J_Phase1_SCG);
figure
plot(xi,f)
title(['SID: ',num2str(SID)])
%%% set two ends of PDF equal to zero
[~,indMode] = max(f);
f(abs(xi-xi(indMode))>80) = 0;
hold on
plot(xi,f)
xlabel('Time [ms]')
%%% Phase 2: find the J-wave with highest PDF
[J_Phase2_SCG, flag_SCG] = FindJ_PDF(J_PEP_SCG,ftSCG_HF.J_A,xi,f);
figure(1)
subplot(6,1,6)
hold on
plot((J_Phase2_SCG.Val(flag_SCG>0,1)+Rpeak(flag_SCG>0))/fs,J_Phase2_SCG.Val(flag_SCG>0,2),'*k')
plot((J_Phase2_SCG.Val(flag_SCG==0,1)+Rpeak(flag_SCG==0))/fs,J_Phase2_SCG.Val(flag_SCG==0,2),'*m')
xticks(RpkKeepInt)
grid on

%% Extract H, I, K, L waves in SCG_HF (chest BCG)------------
ftSCG_HF = FindRestFeat(ftSCG_HF, J_Phase2_SCG.Val(:,1:2), SCGMA_HF, 'SCG',param,Rpeak,fs,1,length(Rpeak));

% Plot extracted H, I, K, L waves in SCG_HF (chest BCG)------------
figure(1)
subplot(6,1,6)
plot(ftSCG_HF.H(:,2)/fs,ftSCG_HF.H(:,1),'om')
plot(ftSCG_HF.I(:,2)/fs,ftSCG_HF.I(:,1),'og')
plot(ftSCG_HF.K(:,2)/fs,ftSCG_HF.K(:,1),'or')
plot(ftSCG_HF.L(:,2)/fs,ftSCG_HF.L(:,1),'ob')
legend('SCG_F_i_l_t','SCG_M_A')

% SCG anterior-posterior feature extraction-----------
ft = FindSCGFeatures(ft,SCGMA,'SCG',param,Rpeak,fs,1,length(Rpeak));
% obj = [];
% obj = scgAortic(obj, varargin); % GT function to extract AO/AC

% Plot anterior posterior SCG features---------
figure(1)
subplot(6,1,5)
plot(ft.AO(:,2)/fs,ft.AO(:,1),'ok')
plot(ft.AC(:,2)/fs,ft.AC(:,1),'*r')

% Calculate and plot BCG J-peak time compared to SCG_HF J-peak and SCG AO-----------
PEP_BCG_SCG = J_Phase2_SCG.Val(and(flag_SCG>0,flag>0),1)-J_Phase2.Val(and(flag_SCG>0,flag>0),1);
PEP_BCG_AO = ft.AO(flag>0,2)-Rpeak(flag>0)-J_Phase2.Val(flag>0,1);
PEP_SCG_AO = ft.AO(and(flag_SCG>0,flag>0),2)-Rpeak(and(flag_SCG>0,flag>0))-J_Phase2_SCG.Val(and(flag_SCG>0,flag>0),1);
figure 
subplot(3,1,1)
plot(PEP_BCG_SCG,'o')
ylabel('J_S_C_G_(_H_F_) - J_B_C_G [ms]')
subplot(3,1,2)
plot(PEP_BCG_AO,'o')
ylabel('AO - J_B_C_G [ms]')
subplot(3,1,3)
plot(PEP_SCG_AO,'o')
ylabel('AO - J_S_C_G_(_H_F_) [ms]')
xlabel('beat #')

% Calculate PTT-------------------------------------
PTT = -J_Phase2.Val(:,1)+ft.tanPPGFt(:,2)-Rpeak(1:end-1); % PTT as the time interval b/w BCG J-peak and PPG foot
PTT_SCG_HF = -J_Phase2_SCG.Val(:,1)+ft.tanPPGFt(:,2)-Rpeak(1:end-1); % PTT as the time interval b/w SCG_HF J-peak and PPG foot
PTT_SCG = -ft.AO(:,2)+ft.tanPPGFt(:,2); % PTT as the time interval b/w anterior posterior SCG AO and PPG foot

%% Find Start Beat # of Each Intervention--------
diff = startIdxCut'-RpkKeep;
[~,startBeat] = min(abs(diff));

% Plot PTT vs. BP Relationship------------------
beatVec = 1:length(Rpeak)-1;
figure
ax3(1) = subplot(2,2,1);
plot(beatVec(~isnan(PTT)),PTT(~isnan(PTT)),'.k')
currAxis = gca;
textloc = 0.9;
PlotInterVLineBeat(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,1,startBeat,length(RpkKeep))
ylabel('PTT [ms]')
ax4(1) = subplot(2,2,2);
plot(beatVec(~isnan(PTT_SCG_HF)),PTT_SCG_HF(~isnan(PTT_SCG_HF)),'.r')
currAxis = gca;
textloc = 0.85;
PlotInterVLineBeat(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,1,startBeat,length(RpkKeep))
ylabel('PTT SCG_H_F [ms]')
ax3(2) = subplot(2,2,3);
plot(beatVec(~isnan(PTT)),ft.tanDP(~isnan(PTT),1),'.k')
currAxis = gca;
textloc = 0.95;
PlotInterVLineBeat(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,1,startBeat,length(RpkKeep))
ylabel('BP [mmHg]')
xlabel('Beat #')
ax4(2) = subplot(2,2,4);
plot(beatVec(~isnan(PTT_SCG_HF)),ft.tanDP(~isnan(PTT_SCG_HF),1),'.r')
currAxis = gca;
textloc = 0.95;
PlotInterVLineBeat(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,1,startBeat,length(RpkKeep))
ylabel('BP [mmHg]')
xlabel('Beat #')
linkaxes(ax3,'x')
linkaxes(ax4,'x')

% Smooth Extracted Features---------------------
Ind = 1:length(ft.tanDP(:,1));
DP_SM = nan(length(ft.tanDP(:,1)),1);  % Preallocate DP matrix
DP_SM(Ind(flag>0)) = ft.tanDP(flag>0,1);  % Fill DP matrix with values from flag>0 beats
SP_SM = nan(length(ft.tanSP(:,1)),1);  % Preallocate SP matrix
SP_SM(Ind(flag>0)) = ft.tanSP(flag>0,1);  % Fill SP matrix with values from flag>0 beats
PTT_SM = nan(length(Ind),1);  % Preallocate PTT matrix
PTT_SM(Ind(flag>0)) = PTT(flag>0,1);  % Fill PTT matrix with values from flag>0 beats
PTT_SCG_HF_SM = nan(length(Ind),1);  % Preallocate PTT_SCG head to foot matrix
PTT_SCG_HF_SM(Ind(and(flag>0,flag_SCG>0))) = PTT_SCG_HF(and(flag>0,flag_SCG>0),1);  % Fill PTT_SCG head to foot matrix with values from flag>0 & flag_SCG>0 beats
PTT_SCG_SM = nan(length(Ind),1);  % Preallocate PTT_SCG matrix
PTT_SCG_SM(Ind(flag>0)) = PTT_SCG(flag>0,1);  % Fill PTT_SCG matrix with values from flag>0 beats
% Omit outlier PTTs with threshold
PTT_SM(PTT_SM<mean(PTT_SM,'omitnan')-std(PTT_SM,'omitnan')) = nan;
PTT_SCG_HF_SM(PTT_SCG_HF_SM<mean(PTT_SCG_HF_SM,'omitnan')-std(PTT_SCG_HF_SM,'omitnan')) = nan;
PTT_SCG_SM(PTT_SCG_SM<mean(PTT_SCG_SM,'omitnan')-std(PTT_SCG_SM,'omitnan')) = nan;
% Hampel and MovMean
HW = 20;
MAW = 20;
DP_SM = hampel(DP_SM,HW);
DP_SM = movmean(DP_SM,MAW);
SP_SM = hampel(SP_SM,HW);
SP_SM = movmean(SP_SM,MAW);
PTT_SM = hampel(PTT_SM,HW);
PTT_SM = movmean(PTT_SM,MAW);
PTT_SCG_HF_SM = hampel(PTT_SCG_HF_SM,HW);
PTT_SCG_HF_SM = movmean(PTT_SCG_HF_SM,MAW);
PTT_SCG_SM = hampel(PTT_SCG_SM,HW);
PTT_SCG_SM = movmean(PTT_SCG_SM,MAW);

% BP Binning-----------------------------------------
[dpBin, temp] = BPBin(DP_SM,PTT_SM);
pttBin(1,:) = temp;
pttBinLable{1} = {'DP PTT_BCG'};
[~, temp] = BPBin(DP_SM,PTT_SCG_HF_SM);
pttBin(2,:) = temp;
pttBinLable{2} = {'DP PTT_SCG_HF'};
[~, temp] = BPBin(DP_SM,PTT_SCG_SM);
pttBin(3,:) = temp;
pttBinLable{3} = {'DP PTT_SCG'};
[spBin, temp] = BPBin(SP_SM,PTT_SM);
pttBin_SP(1,:) = temp;
pttBin_SPLable{1} = {'SP PTT_BCG'};
[~, temp] = BPBin(SP_SM,PTT_SCG_HF_SM);
pttBin_SP(2,:) = temp;
pttBin_SPLable{2} = {'SP PTT_SCG_HF'};
[~, temp] = BPBin(SP_SM,PTT_SCG_SM);
pttBin_SP(3,:) = temp;
pttBin_SPLable{3} = {'SP PTT_SCG'};
%%% Regression Bin
TableCorr = [];
x = pttBin(1,:);  
y = dpBin;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{1,1} = {'Bin: DP vs. PTT_BCG'};
x = pttBin(2,:);  
y = dpBin;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{2,1} = {'Bin: DP vs. PTT_SCG_HF'};
x = pttBin(3,:);  
y = dpBin;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{3,1} = {'Bin: DP vs. PTT_SCG'};
x = pttBin_SP(1,:);  
y = spBin;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{4,1} = {'Bin: SP vs. PTT_BCG'};
x = pttBin_SP(2,:);  
y = spBin;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{5,1} = {'Bin: SP vs. PTT_SCG_HF'};
x = pttBin_SP(3,:);  
y = spBin;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{6,1} = {'Bin: SP vs. PTT_SCG'};

% Regression All points----------------------------
x = PTT_SM;  
y = DP_SM;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{7,1} = {'All: DP vs. PTT_BCG'};
x = PTT_SCG_HF_SM;  
y = DP_SM;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{8,1} = {'All: DP vs. PTT_SCG_HF'};
x = PTT_SCG_SM;  
y = DP_SM;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{9,1} = {'All: DP vs. PTT_SCG'};
x = PTT_SM;  
y = SP_SM;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{10,1} = {'All: SP vs. PTT_BCG'};
x = PTT_SCG_HF_SM;  
y = SP_SM;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{11,1} = {'All: SP vs. PTT_SCG_HF'};
x = PTT_SCG_SM;  
y = SP_SM;  
mdl = fitlm(x,y,'linear');
R = sqrt(mdl.Rsquared.Ordinary);
slope = mdl.Coefficients{2,1};
intrcpt = mdl.Coefficients{1,1};
TableCorr = [TableCorr;R slope intrcpt];
TableCorrLabel{12,1} = {'All: SP vs. PTT_SCG'};

figure
subplot(2,3,1)
plot(pttBin(1,:),dpBin,'o')
hold on
plot(pttBin(1,:),TableCorr(1,2)*pttBin(1,:)+TableCorr(1,3))
ylabel('DP')
title('BCG')
subplot(2,3,2)
plot(pttBin(2,:),dpBin,'o')
hold on
plot(pttBin(2,:),TableCorr(2,2)*pttBin(2,:)+TableCorr(2,3))
title('SCG_H_F')
subplot(2,3,3)
plot(pttBin(3,:),dpBin,'o')
hold on
plot(pttBin(3,:),TableCorr(3,2)*pttBin(3,:)+TableCorr(3,3))
title('SCG')
subplot(2,3,4)
plot(pttBin_SP(1,:),spBin,'o')
hold on
plot(pttBin_SP(1,:),TableCorr(4,2)*pttBin_SP(1,:)+TableCorr(4,3))
ylabel('SP')
xlabel('PTT [ms]')
subplot(2,3,5)
plot(pttBin_SP(2,:),spBin,'o')
hold on
plot(pttBin_SP(2,:),TableCorr(5,2)*pttBin_SP(2,:)+TableCorr(5,3))
xlabel('PTT [ms]')
subplot(2,3,6)
plot(pttBin_SP(3,:),spBin,'o')
hold on
plot(pttBin_SP(3,:),TableCorr(6,2)*pttBin_SP(3,:)+TableCorr(6,3))
xlabel('PTT [ms]')

if SID == 1
    axLim = [50 110 170;180 140 230];
elseif SID == 2
    axLim = [60 80 170;140 130 280];
elseif SID == 3
    axLim = [70 100 170;170 140 250];
elseif SID == 4
    axLim = [70 100 170;170 140 250];
elseif SID == 5
    axLim = [70 100 170;170 140 250];
end

figure
ax2(1) = subplot(5,1,1);
plot(PTT_SM,'.')
ylabel('PTT [ms] (BCG)')
axis([-inf inf axLim(1) axLim(2)])
currAxis = gca;
textloc = 0.8;
PlotInterVLineBeat(currAxis.YLim(2),timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,1,startBeat,length(RpkKeep))
ax2(2) = subplot(5,1,2);
plot(PTT_SCG_HF_SM,'.')
ylabel('PTT [ms] (SCG_H_F)')
axis([-inf inf axLim(3) axLim(4)])
ax2(3) = subplot(5,1,3);
plot(PTT_SCG_SM,'.')
ylabel('PTT [ms] (SCG)')
axis([-inf inf axLim(5) axLim(6)])
ax2(4) = subplot(5,1,4);
plot(DP_SM,'.r')
ylabel('DP [mmHg]')
ax2(5) = subplot(5,1,5);
plot(SP_SM,'.r')
ylabel('SP [mmHg]')
linkaxes(ax2,'x')

figure
subplot(2,3,1)
plot(PTT_SM,DP_SM,'.')
hold on
plot(PTT_SM(~isnan(DP_SM)),TableCorr(7,2)*PTT_SM(~isnan(DP_SM))+TableCorr(7,3))
ylabel('DP [mmHg]')
title('BCG')
plot(pttBin(1,:),dpBin,'*k')
subplot(2,3,2)
plot(PTT_SCG_HF_SM,DP_SM,'.')
hold on
plot(PTT_SCG_HF_SM(~isnan(DP_SM)),TableCorr(8,2)*PTT_SCG_HF_SM(~isnan(DP_SM))+TableCorr(8,3))
plot(pttBin(2,:),dpBin,'*k')
title('SCG_H_F')
subplot(2,3,3)
plot(PTT_SCG_SM,DP_SM,'.')
hold on
plot(PTT_SCG_SM(~isnan(DP_SM)),TableCorr(9,2)*PTT_SCG_SM(~isnan(DP_SM))+TableCorr(9,3))
plot(pttBin(3,:),dpBin,'*k')
title('SCG')
subplot(2,3,4)
plot(PTT_SM,SP_SM,'.r')
hold on
plot(PTT_SM(~isnan(SP_SM)),TableCorr(10,2)*PTT_SM(~isnan(SP_SM))+TableCorr(10,3))
plot(pttBin_SP(1,:),spBin,'*k')
ylabel('SP [mmHg]')
xlabel('PTT [ms]')
subplot(2,3,5)
plot(PTT_SCG_HF_SM,SP_SM,'.r')
hold on
plot(PTT_SCG_HF_SM(~isnan(SP_SM)),TableCorr(11,2)*PTT_SCG_HF_SM(~isnan(SP_SM))+TableCorr(11,3))
plot(pttBin_SP(2,:),spBin,'*k')
xlabel('PTT [ms]')
subplot(2,3,6)
plot(PTT_SCG_SM,SP_SM,'.r')
hold on
plot(PTT_SCG_SM(~isnan(SP_SM)),TableCorr(12,2)*PTT_SCG_SM(~isnan(SP_SM))+TableCorr(12,3))
plot(pttBin_SP(3,:),spBin,'*k')
xlabel('PTT [ms]')

% Tables-----------------------------------------------
TablePercentFeature = [sum(flag>0)/length(Rpeak)*100, sum(flag_SCG>0)/length(Rpeak)*100];  % percentage BCG, SCG_HF beats with extracted J-peak
% Mean and SD Values
TableValLabels = {'DP';'SP';'PTT_BCG';'PTT_SCG_HF';'PTT_SCG';'PEP_BCG-PEP_SCG_HF';'PEP_BCG-AO';'PEP_SCG_HF-AO'};
TableVal(1,:) = [mean(DP_SM,'omitnan'), std(DP_SM,'omitnan'), min(DP_SM), max(DP_SM)];
TableVal(2,:) = [mean(SP_SM,'omitnan'), std(SP_SM,'omitnan'), min(SP_SM), max(SP_SM)];
TableVal(3,:) = [mean(PTT_SM,'omitnan'), std(PTT_SM,'omitnan'), min(PTT_SM), max(PTT_SM)];
TableVal(4,:) = [mean(PTT_SCG_HF_SM,'omitnan'), std(PTT_SCG_HF_SM,'omitnan'), min(PTT_SCG_HF_SM), max(PTT_SCG_HF_SM)];
TableVal(5,:) = [mean(PTT_SCG_SM,'omitnan'), std(PTT_SCG_SM,'omitnan'), min(PTT_SCG_SM), max(PTT_SCG_SM)];
TableVal(6,:) = [mean(PEP_BCG_SCG), std(PEP_BCG_SCG), min(PEP_BCG_SCG), max(PEP_BCG_SCG)];
TableVal(7,:) = [mean(PEP_BCG_AO), std(PEP_BCG_AO), min(PEP_BCG_AO), max(PEP_BCG_AO)];
TableVal(8,:) = [mean(PEP_SCG_AO), std(PEP_SCG_AO), min(PEP_SCG_AO), max(PEP_SCG_AO)];
