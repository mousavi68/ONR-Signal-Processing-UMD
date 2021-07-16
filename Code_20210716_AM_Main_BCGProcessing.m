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

% Add Path to library if needed-------------------
currentfolder = pwd;

% Load Excel file: Intervention Info---------------
[~,~,subjectInterventionXLS] = xlsread('subjectIntervention'); % Excel file containing interventions timings and info on whether or not include each intervention

% Load Data-------------------------------------------
load('Weidi_full_protocol_05202021'); SID = 11; fileIndex = 1; % Change SID and fileIndex for each new subject. fileIndex should be set equl to the line # in the excel file containing current subject info subtracted by 1

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
warning('off') % turn warnings off

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
RpkKeep = Rpeak; % ECG R-peak in sample #
RpkKeepInt = Rpeak/fs; % ECG R-peak timing in [s]

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

% BP Processing and feature extraction---------
[u,v] = butter(2,20/(fs/2));
BP = filtfilt(u,v,BP);
ft = []; % Variable storing the features
ft = FindBPFeatures(ft, BP,'BP',param.BP,Rpeak,fs,1,length(Rpeak)); % updating ft with BP features

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
ft = FindBPFeatures(ft, PPG,'PPG',param.PPG,Rpeak,fs,1,length(Rpeak));

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

% SCG Preprocessing (?)------------------------------
[u,v]=butter(2,[0.5,10]/(fs/2));
% SCG = filtfilt(u,v,SCG);
SCGBut_HF = filtfilt(u,v,SCG_HF);
% GT Filter for SCG and BCG
Hd = createBPF(0.5,40,fs,'fpass1',1,'fpass2', 20, 'kaiser', 'order', 20);
BCGKsr_W = filtfilt(Hd.Numerator,1,BCG_W);
SCGFilt = filtfilt(Hd.Numerator, 1, SCG);
% SCGFilt_HF = filtfilt(Hd.Numerator, 1, SCG_HF);
SCGFilt_HF = SCGBut_HF;

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
%%% reference for this section (PDF based J-peak detection): S. Shin et. al., "A Unified Approach to
%%% Wearable Ballistocardiogram Gating and Wave Localization", 2020
%%% Phase 1: find j-wave candidates, each beat can have 10 candidates at the most
ft = FindFeatures(ft, BCGMA_W, 'BCG',param.BCG_A,Rpeak,fs,1,length(Rpeak)); % Find J-wave candidates 
J_PEP = ft.J_T-ones(10,1)*transpose(Rpeak(1:end-1));
[~,ind] = max(ft.J_A);
for ii = 1:length(ft.J_A)
J_Phase1(ii) = J_PEP(ind(ii),ii);
end
[f,xi] = ksdensity(J_Phase1); % matlab's built-in function for probability density estimate
%%% set two ends that are further than 80 [ms] from the max of PDF equal to zero
[~,indMode] = max(f);
f(abs(xi-xi(indMode))>80) = 0;  
%%% Phase 2: find the J-wave with highest PDF
[J_Phase2, flag] = FindJ_PDF(J_PEP,ft.J_A,xi,f);

% Extract H, I, K, L waves in wrist BCG------------
ft = FindRestFeat(ft, J_Phase2.Val(:,1:2), BCGMA_W, 'BCG_A',param,Rpeak,fs,1,length(Rpeak));

% Plot extracted BCG Features-------------------
figure(1)
subplot(6,1,4)
hold on
plot((J_Phase2.Val(flag>0,1)+Rpeak(flag>0))/fs,J_Phase2.Val(flag>0,2),'*k')
xticks(RpkKeepInt)
plot(ft.H(:,2)/fs,ft.H(:,1),'om')
plot(ft.I(:,2)/fs,ft.I(:,1),'og')
plot(ft.K(:,2)/fs,ft.K(:,1),'or')
plot(ft.L(:,2)/fs,ft.L(:,1),'ob')
grid on

% Chest SCG Head-to-Foot (Chest BCG) J-wave Extraction--------
%%% Phase 1: find j-wave candidates, each beat can have 10 candidates max
ftSCG_HF = [];
[ftSCG_HF,plotBeat] = FindFeatures(ftSCG_HF, SCGMA_HF, 'BCG',param.SCG,Rpeak,fs,1,length(Rpeak));
% figure(1)
% subplot(6,1,6)
% plot(ftSCG_HF.J_T/fs,ftSCG_HF.J_A,'o')
J_PEP_SCG = ftSCG_HF.J_T-ones(10,1)*transpose(Rpeak(1:end-1));
[~,ind] = max(ftSCG_HF.J_A);
for ii = 1:length(ftSCG_HF.J_A)
J_Phase1_SCG(ii) = J_PEP_SCG(ind(ii),ii);
end
[f,xi] = ksdensity(J_Phase1_SCG);
%%% set two ends of PDF equal to zero
[~,indMode] = max(f);
f(abs(xi-xi(indMode))>80) = 0;
%%% Phase 2: find the J-wave with highest PDF
[J_Phase2_SCG, flag_SCG] = FindJ_PDF(J_PEP_SCG,ftSCG_HF.J_A,xi,f);
figure(1)
subplot(6,1,6)
hold on
plot((J_Phase2_SCG.Val(flag_SCG>0,1)+Rpeak(flag_SCG>0))/fs,J_Phase2_SCG.Val(flag_SCG>0,2),'*k')
% plot((J_Phase2_SCG.Val(flag_SCG==0,1)+Rpeak(flag_SCG==0))/fs,J_Phase2_SCG.Val(flag_SCG==0,2),'*m')
xticks(RpkKeepInt)
grid on

% Extract H, I, K, L waves in SCG_HF (chest BCG)------------
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

% Plot anterior posterior SCG features---------
figure(1)
subplot(6,1,5)
plot(ft.AO(:,2)/fs,ft.AO(:,1),'ok')
plot(ft.AC(:,2)/fs,ft.AC(:,1),'*r')

% Calculate BCG J-peak time compared to SCG_HF J-peak and SCG AO-----------
PEP_BCG_SCG = J_Phase2_SCG.Val(and(flag_SCG>0,flag>0),1)-J_Phase2.Val(and(flag_SCG>0,flag>0),1);
PEP_BCG_AO = ft.AO(flag>0,2)-Rpeak(flag>0)-J_Phase2.Val(flag>0,1);
PEP_SCG_AO = ft.AO(and(flag_SCG>0,flag>0),2)-Rpeak(and(flag_SCG>0,flag>0))-J_Phase2_SCG.Val(and(flag_SCG>0,flag>0),1);

% Calculate PTT-------------------------------------
PTT = -J_Phase2.Val(:,1)+ft.tanPPGFt(:,2)-Rpeak(1:end-1); % PTT as the time interval b/w BCG J-peak and PPG foot
PTT_SCG_HF = -J_Phase2_SCG.Val(:,1)+ft.tanPPGFt(:,2)-Rpeak(1:end-1); % PTT as the time interval b/w SCG_HF J-peak and PPG foot
PTT_SCG = -ft.AO(:,2)+ft.tanPPGFt(:,2); % PTT as the time interval b/w anterior posterior SCG AO and PPG foot

% Plot PTT vs. BP Relationship----------------------
figure
ax3(1) = subplot(2,2,1);
plot(PTT(flag>0),'.k')
ylabel('PTT [ms]')
ax4(1) = subplot(2,2,2);
plot(PTT_SCG_HF(flag_SCG>0),'.r')
ylabel('PTT SCG_H_F [ms]')
ax3(2) = subplot(2,2,3);
plot(ft.tanDP(flag>0,1),'.k')
ylabel('BP [mmHg]')
xlabel('Beat #')
ax4(2) = subplot(2,2,4);
plot(ft.tanDP(flag_SCG>0,1),'.r')
ylabel('BP [mmHg]')
xlabel('Beat #')
linkaxes(ax3,'x')
linkaxes(ax4,'x')
