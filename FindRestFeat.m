function [ft,plotBeat] = FindRestFeat(ft, J, signal, signalName,param,peak,fs,stBt,ndBt)
% Description:-----------------------------------------
% This function extracts H, I, K, L waves in the BCG/SCG head-to-foot
% The function uses J-peaks in each beat as a starting point. It goes back
% from the J-peak and searches the area specified by ISt to location of
% J-peak to extract I-peak. If I-peak exists in the current beat, the
% algorithm searches the area specified by HSt to location of I-wave to
% find H-wave that is larger than a threshold value specified by
% HMinPeakHeight. Similar procedure is used for extracting K and L wave
% with one difference. To extract K-wave the algorithm starts from J-wave
% location and goes forward and seraches the area upto KEnd for a candidate
% that is larger than KMinPeakHeight. If K-wave exists the algorithm
% seraches the area after K-wave upto LEnd for a candidate that is larger
% than LMinPeakHeight. In case J doesn't exist in the current beat all
% other features will be set to Nan. In case, I doesn't exists in the
% current beat I and H are set to Nan, and so on.

% Inputs:------------------------------------------------
% ft: variable storing features
% J: contains the location and amplitude of the extracted J-waves
% signal: filtered BCG or SCG_HF
% signalName: set to 'BCG_A' for BCG and set to 'SCG' for SCG_HF
% param: timing and amplitude values used in findpeaks function to extract the features
% peak: ECG R-peaks
% fs: sampling frequency
% stBt: start beat # to extract features
% ndBt: last beat # 

% Outputs:---------------------------------------------
% ft: updated variable storing features
% plotBeat: ECG R-peaks not including NaN

% Body of the Code:---------------------------------
% beatNo = length(peak)-1;
plotBeat = ~isnan(peak(1:end-1));
num_I = 0; num_K = 0; num_H = 0; num_L = 0;
for ii = stBt:ndBt-1
        if ~isnan(peak(ii)) && ~isnan(peak(ii+1)) && ~isempty(J(ii,1)) && ~isnan(J(ii,1))
                % I wave
                if ~isempty(findpeaks(-flip(signal(peak(ii)+param.(signalName).ISt:J(ii,1)+peak(ii))),'minpeakheight',param.(signalName).IMinPeakHeight,'NPeaks',1))
                    [ft.I(ii,1),tempLoc] = findpeaks(-flip(signal(peak(ii)+param.(signalName).ISt:J(ii,1)+peak(ii))),'minpeakheight',param.(signalName).IMinPeakHeight,'NPeaks',1);
                    ft.I(ii,1) = -ft.I(ii,1);
                    ft.I(ii,2) = length(signal(peak(ii)+param.(signalName).ISt:J(ii,1)+peak(ii)))-tempLoc+peak(ii)+param.(signalName).ISt;
                    num_I = num_I+1;
                    Ind.I(num_I) = ii;
                    % H wave
                    if ~isempty(findpeaks(flip(signal(peak(ii)+param.(signalName).HSt:ft.I(ii,2))),'minpeakheight',param.(signalName).HMinPeakHeight,'NPeaks',1))
                        [ft.H(ii,1),tempLoc] = findpeaks(flip(signal(peak(ii)+param.(signalName).HSt:ft.I(ii,2))),'minpeakheight',param.(signalName).HMinPeakHeight,'NPeaks',1);
                        ft.H(ii,1) = ft.H(ii,1);
                        ft.H(ii,2) = length(signal(peak(ii)+param.(signalName).HSt:ft.I(ii,2)))-tempLoc+peak(ii)+param.(signalName).HSt;
                        num_H = num_H+1;
                        Ind.H(num_H) = ii;
                    else
                        ft.H(ii,1:2) = [nan peak(ii)];
                    end
                else
                    ft.I(ii,1:2) = [nan peak(ii)];
                    ft.H(ii,1:2) = [nan peak(ii)];
                end
                % K wave
                if ~isempty(findpeaks(-signal(J(ii,1)+peak(ii):J(ii,1)+param.(signalName).KEnd+peak(ii)),'minpeakheight',param.(signalName).KMinPeakHeight,'NPeaks',1))
                    [ft.K(ii,1),tempLoc] = findpeaks(-signal(J(ii,1)+peak(ii):J(ii,1)+param.(signalName).KEnd+peak(ii)),'minpeakheight',param.(signalName).KMinPeakHeight,'NPeaks',1);
                    ft.K(ii,1) = -ft.K(ii,1);
                    ft.K(ii,2) = tempLoc+J(ii,1)+peak(ii);
                    num_K = num_K+1;
                    Ind.K(num_K) = ii;
                    % L wave
                    if ~isempty(findpeaks(signal(ft.K(ii,2):ft.K(ii,2)+param.(signalName).LEnd),'minpeakheight',param.(signalName).LMinPeakHeight,'NPeaks',1))
                        [ft.L(ii,1),tempLoc] = findpeaks(signal(ft.K(ii,2):ft.K(ii,2)+param.(signalName).LEnd),'minpeakheight',param.(signalName).LMinPeakHeight,'NPeaks',1);
                        ft.L(ii,1) = ft.L(ii,1);
                        ft.L(ii,2) = tempLoc+ft.K(ii,2);
                        num_L = num_L+1;
                        Ind.L(num_L) = ii;
                    else
                         ft.L(ii,1:2) = [nan peak(ii)];
                    end
                else
                    ft.K(ii,1:2) = [nan peak(ii)];
                    ft.L(ii,1:2) = [nan peak(ii)];
                end
        else
            plotBeat(ii) = 0;
            ft.I(ii,1:2) = [nan peak(ii)];
            ft.K(ii,1:2) = [nan peak(ii)];
            ft.H(ii,1:2) = [nan peak(ii)];
            ft.L(ii,1:2) = [nan peak(ii)];
        end
end
    if exist('Ind')
        if ~isfield(Ind,'I')
            Ind.I = [];
        end
    end
end
