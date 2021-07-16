function [flagBeat,flagVal,flagInd,peakExc] = AssessBeats(accSig,peaks,mult)
% Description:-----------------------------------------
% This function excludes BCG beats with large motions

% Input:-------------------------------------------------
% accSig: cardiomechanical signal under investigation (BCG or SCG)
% peaks: ECG R-peaks
% mult: controls the threshold for large motion based beat exclusion 

% Output:-----------------------------------------------
% flagBeat: 1 or 0. 1 means all/some portion of the signal in the current
% beat is larger than the threshold and therefore the entire beat should be
% discarded
% flagVal: absolut value of the threshold
% flagInd: sample # which will be excluded due to motion artifact
% peakExc: remaining ECG R-peaks after implementing the motion artifact
% exclusion criteria

% Body of the Code:--------------------------------
peakExc = peaks;
meanSig = mean(abs(accSig));
medSig = median(abs(accSig));
flag = zeros(1,length(accSig));
flagInd =[];
flagVal = mult*medSig;
flag(abs(accSig)>flagVal) = 1;
% figure
% plot(accSig)
% hold on
% plot(flagVal*ones(1,length(accSig)))
% flag(accSig>5*meanSig) = 1;
for ii = 1:length(peaks)-1
    if sum(flag(peaks(ii):peaks(ii+1)))>0
        flagBeat(ii,1) = 1;
        peakExc(ii) = nan;
    else
        flagBeat(ii,1) = 0;
    end
end
% plot(flagInd,accSig(flagInd),'.')
% flagBtBef = flagBeat(2:end-3)+flagBeat(1:end-4);
% flagBtAft = flagBeat(4:end-1)+flagBeat(5:end);
% flagBtBef = [0,0,flagBtBef,0,0];
% flagBtAft = [0,0,flagBtAft,0,0];
flagBtBef = flagBeat(1:end-2);
flagBtAft = flagBeat(3:end);
flagBtBef = [0;flagBtBef;0];
flagBtAft = [0;flagBtAft;0];
flagBeat((flagBtBef>0)&(flagBtAft>0)) = 1;
peakExc(flagBeat==1) = nan;

for ii = 1:length(peaks)-1
    if flagBeat(ii) == 1
         flagInd = [flagInd,(peaks(ii):peaks(ii+1)-1)];  % sample # with large motion
    end
end
end