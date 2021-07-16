function [ft,plotBeat] = FindBPFeatures(ft, signal, signalName,param,peak,fs,stBt,ndBt)
% Description:-----------------------------------------
% This function extracts Systolic and diastolic pressures (as well as PPG min/max) via two
% approaches:1. min/max  of BP in each beat 2. intersecting tangent method

% Inputs:------------------------------------------------
% ft: variable containing features
% signal: filtered BP or PPG signal
% signal name: set to 'BP' or 'PPG'
% param: values of parameters used in function findpeaks
% fs: sampling frequency
% stBt: starting beat #
% ndBt: last beat number

% Outputs: --------------------------------------------
% ft: updated variable storing features
% plotBeat: ECG R-peaks not including NaN

% Body of the Code: --------------------------------
% beatNo = length(peak)-1;
plotBeat = ~isnan(peak(1:end-1));
if strcmp(signalName,'BP')
    num_BP = 0;
    for ii = stBt:ndBt-1
        if ii == ndBt-1
            ii;
        end
        if ~isnan(peak(ii)) && ~isnan(peak(ii+1)) && ~isempty(findpeaks(signal(peak(ii)+param.SPSt:peak(ii)+param.SPEnd),'minpeakheight',param.SPMinPeakHeight))
            % Find SP/DP based on min/max approach
            [ft.SP(ii,1),tempLoc] = findpeaks(signal(peak(ii)+param.SPSt:peak(ii)+param.SPEnd),'minpeakheight',param.SPMinPeakHeight,'NPeaks',1);
            ft.SP(ii,2) = tempLoc+peak(ii)+param.SPSt;
            if ~isempty(findpeaks(-signal(peak(ii):ft.SP(ii,2)),'minpeakheight',param.DPMinPeakHeight,'NPeaks',1))
                [ft.DP(ii,1),tempLoc] = findpeaks(-signal(peak(ii):ft.SP(ii,2)),'minpeakheight',param.DPMinPeakHeight,'NPeaks',1);
                ft.DP(ii,2) = tempLoc+peak(ii);
                ft.DP(ii,1) = -ft.DP(ii,1);
                % Find SP/DP based on Intersecting Tangent Method
                [ft.tanDP(ii,2),ft.tanSP(ii,2)] = IntTang(signal,ft.DP(ii,2),ft.SP(ii,2));
                ft.tanDP(ii,1) = ft.DP(ii,1);
                ft.tanSP(ii,1) = ft.SP(ii,1);
            else
                ft.DP(ii,1:2) = [nan peak(ii)];
                ft.tanDP(ii,1:2) = [nan peak(ii)];
                ft.tanSP(ii,1:2) = [nan peak(ii)];
            end
            num_BP = num_BP+1;
            Ind(num_BP) = ii;
        else
            ft.SP(ii,1:2) = [nan peak(ii)];
            ft.DP(ii,1:2) = [nan peak(ii)];
            ft.tanDP(ii,1:2) = [nan peak(ii)];
            ft.tanSP(ii,1:2) = [nan peak(ii)];
        end
    end
elseif strcmp(signalName,'PPG')
    for ii = stBt:ndBt-1
        if ~isnan(peak(ii)) && ~isnan(peak(ii+1)) && ~isempty(findpeaks(signal(peak(ii)+param.PPGPkSt:peak(ii)+param.PPGPkEnd),'minpeakheight',param.PPGMinPeakHeight))
            % Find SP/DP based on min/max approach
            [ft.PPGPk(ii,1),tempLoc] = findpeaks(signal(peak(ii)+param.PPGPkSt:peak(ii)+param.PPGPkEnd),'minpeakheight',param.PPGMinPeakHeight,'NPeaks',1);
            ft.PPGPk(ii,2) = tempLoc+peak(ii)+param.PPGPkSt;
            if ~isempty(findpeaks(-signal(peak(ii):ft.PPGPk(ii,2)),'minpeakheight',param.PPGFtMinPeakHeight,'NPeaks',1))
                [ft.PPGFt(ii,1),tempLoc] = findpeaks(-signal(peak(ii):ft.PPGPk(ii,2)),'minpeakheight',param.PPGFtMinPeakHeight,'NPeaks',1);
                ft.PPGFt(ii,2) = tempLoc+peak(ii);
                ft.PPGFt(ii,1) = -ft.PPGFt(ii,1);
                % Find Peak/Foot based on Intersecting Tangent Method
                [ft.tanPPGFt(ii,2),ft.tanPPGPk(ii,2)] = IntTang(signal,ft.PPGFt(ii,2),ft.PPGPk(ii,2));
                ft.tanPPGFt(ii,1) = ft.PPGFt(ii,1);
                ft.tanPPGPk(ii,1) = ft.PPGPk(ii,1);
            else
                ft.PPGFt(ii,1:2) = [nan peak(ii)];
                ft.tanPPGFt(ii,1:2) = [nan peak(ii)];
                ft.tanPPGPk(ii,1:2) = [nan peak(ii)];
            end
        else
            ft.PPGPk(ii,1:2) = [nan peak(ii)];
            ft.PPGFt(ii,1:2) = [nan peak(ii)];
            ft.tanPPGFt(ii,1:2) = [nan peak(ii)];
            ft.tanPPGPk(ii,1:2) = [nan peak(ii)];
        end
    end
end
end