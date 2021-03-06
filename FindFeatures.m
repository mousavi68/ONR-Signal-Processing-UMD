function [ft,plotBeat] = FindFeatures(ft, signal, signalName,param,peak,fs,stBt,ndBt)
% Description: ----------------------------------------
% function extracts multiple J-peak candidates in each beat for further
% analysis. Matlab's bult-in function "findpeaks" is used to search for candidates in an
% area specified by parameter values JSt as the starting point and JEnd as
% the end point

% Inputs:-----------------------------------------------
% ft: variable containing features
% signal: BCG
% signalName: set to 'BCG'
% param: defined in the function SetParameters, contains the values used in
% the findpeaks function for extracting J-wave candidates
% peak: ECG R-peaks
% fs: sampling frequency
% stBt: starting beat #
% ndBt: last beat #

% Outputs: --------------------------------------------
% ft: updated variable storing features
% plotBeat: ECG R-peaks not including NaN

% Body of the Code:---------------------------------
% beatNo = length(peak)-1;
plotBeat = ~isnan(peak(1:end-1));
N = 10;
% Preallocating the vectors -  a maximum of 10 J-wave candidates is considered 
ft.J_A= nan(N,1);
ft.J_T= nan(N,1);
if strcmp(signalName,'BCG')
    num_J = 0; num_I = 0; num_K = 0; num_H = 0; num_L = 0;
    for ii = stBt:ndBt-1
        if ~isnan(peak(ii)) && ~isnan(peak(ii+1)) && ~isempty(findpeaks(signal(peak(ii)+param.JSt:peak(ii)+param.JEnd),'minpeakheight',param.JMinPeakHeight,'NPeaks',1,'MinPeakProminence',param.JProm))
            % J wave
                [tempA,tempLoc] = findpeaks(signal(peak(ii)+param.JSt:peak(ii)+param.JEnd),'minpeakheight',param.JMinPeakHeight,'MinPeakProminence',param.JProm);
                ft.J_A(:,ii) = [tempA;nan(N-length(tempA),1)];
                ft.J_T(:,ii) = [tempLoc+param.JSt+peak(ii);nan(N-length(tempA),1)];
                num_J = num_J+1;
                Ind.J(num_J) = ii;
        else
            ft.J_A(1:N,ii) = nan(N,1);
            ft.J_T(1:N,ii) = nan(N,1);
            plotBeat(ii) = 0;
        end
    end
end
end