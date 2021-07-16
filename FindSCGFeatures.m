function ft = FindSCGFeatures(ft, signal, signalName,param,peak,fs,stBt,ndBt)
% Description: ----------------------------------------
% function extracts multiple J-peak candidates in each beat for further
% analysis. Matlab's bult-in function "findpeaks" is used to search for candidates in an
% area specified by parameter values JSt as the starting point and JEnd as
% the end point

% Inputs:-----------------------------------------------
% ft: variable containing features
% signal: SCG_HF
% signalName: set to 'SCG'
% param: defined in the function SetParameters, contains the values used in
% the findpeaks function for extracting J-wave candidates
% peak: ECG R-peaks
% fs: sampling frequency
% stBt: starting beat #
% ndBt: last beat #

% Outputs: --------------------------------------------
% ft: updated variable storing features

% Body of the Code:---------------------------------
% beatNo = length(peak)-1;
    if strcmp(signalName,'SCG')
        num_BP = 0;
        for ii = stBt:ndBt-1
            if ~isnan(peak(ii)) && ~isnan(peak(ii+1)) && ~isempty(findpeaks(signal(peak(ii)+param.AOSt:peak(ii)+param.AOEnd),'minpeakheight',param.AOMinPeakHeight))
                % Find SP/DP based on min/max approach
                [ft.AO(ii,1),tempLoc] = findpeaks(signal(peak(ii)+param.AOSt:peak(ii)+param.AOEnd),'minpeakheight',param.AOMinPeakHeight,'NPeaks',1);
                ft.AO(ii,2) = tempLoc+peak(ii)+param.AOSt;
                if ~isempty(findpeaks(signal(ft.AO(ii,2)+param.ACSt:ft.AO(ii,2)+param.ACEnd),'minpeakheight',param.ACMinPeakHeight))
                    [ft.AC(ii,1),tempLoc] = findpeaks(signal(ft.AO(ii,2)+param.ACSt:ft.AO(ii,2)+param.ACEnd),'minpeakheight',param.ACMinPeakHeight,'NPeaks',1);
                    ft.AC(ii,2) = tempLoc+ft.AO(ii,2)+param.ACSt;
                else
                    ft.AC(ii,1:2) = [nan peak(ii)];
                end
                num_BP = num_BP+1;
                Ind(num_BP) = ii;
            else
                ft.AO(ii,1:2) = [nan peak(ii)];
                ft.AC(ii,1:2) = [nan peak(ii)];
            end
        end
    end
end