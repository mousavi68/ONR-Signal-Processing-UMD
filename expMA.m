function BCGsig=expMA(BCGsig,BCGCopysig,ECGMaxLoc,startIdx,time,fs)
% Description:-----------------------------------------
% 10-beat exponential moving average filter for BCG/SCG signal

% Inputs:-----------------------------------------------
% BCGsig: input signal
% BCGCopysig: copy of the input signal 
% ECGMaxLoc: ECG R-peak sample number
% startIdx: Start Index
% time: time signal
% fs: sampling frequency

% Outputs:---------------------------------------------
% BCGsig: updated signal after applying exponential moving average

% Body of the Code:---------------------------------
M=10;alpha=2/(M+1);
ECGMaxDiff=diff(ECGMaxLoc);
nInt=round(nanmedian(ECGMaxDiff)*1.2);
maxPeakIdx=length(ECGMaxLoc);
while (ECGMaxLoc(maxPeakIdx)-(startIdx-1)+nInt-1>length(BCGsig))
    maxPeakIdx = maxPeakIdx-1;
end
% disp(maxPeakIdx)
if maxPeakIdx<M
    fprintf('Error: Not enough beats for BCG exponential moving averaging!\n')
else
    for peakIdx = M:maxPeakIdx
        tempNum=0;tempDen=0;
        for beatIdx = 1:M
            start = ECGMaxLoc(peakIdx-beatIdx+1)-(startIdx-1);
            if ~isnan(start)
%                 [nInt,start,start+nInt,length(BCGsig)]
                tempBCGsig = BCGCopysig(start:start+nInt-1);
                tempAmp(beatIdx,1) = max(tempBCGsig)-min(tempBCGsig);
                isValid(beatIdx,1) = true;
            else
                tempAmp(beatIdx,1) = NaN;
                isValid(beatIdx,1) = false;
            end
        end
%         tempAmp
%         isValid
        isIncluded = double((tempAmp<=3*nanmedian(tempAmp) & isValid));
        for beatIdx = 1:M
            start = ECGMaxLoc(peakIdx-beatIdx+1)-(startIdx-1);
            if isIncluded(beatIdx,1) == 1
                tempNum=tempNum+(1-alpha)^(beatIdx-1)*BCGCopysig(start:start+nInt-1);
                tempDen=tempDen+(1-alpha)^(beatIdx-1);
            end
        end
        newstart = ECGMaxLoc(peakIdx)-(startIdx-1);
        if ~isnan(newstart)
            BCGsig(newstart:newstart+nInt-1)=tempNum/tempDen;
        end
    end
end
% % Integrate- Azin
% % time = getappdata(mPlots.Analysis.EventFiltered.Parent,'timesig');
% % load mparams
% % fs=mParams.Constant.fs; % sampling frequency
% % ts = 1/fs;
% % t1=0:ts:ts*(length(ECG_T3)-1);
% int1 = cumtrapz(time(startIdx:startIdx+length(BCGsig)-1)',BCGsig);
% int2 = cumtrapz(time(startIdx:startIdx+length(BCGsig)-1)',int1);
% global N
% [u,v]=butter(2,[N,50]/(fs/2));
% % [u,v]=butter(2,[0.5,50]/(fs/2));   % Azin
% BCGsig=filtfilt(u,v,int2)*1e3;
% 
% % setappdata(mPlots.Analysis.(eventTag).Parent,'BCGsig',BCGsig)
% % set(mPlots.Analysis.(eventTag).BCG,'YData',BCGsig(currentWindowStartIdx:currentWindowEndIdx));
fprintf(['EXP MA BCG has been calculated!\n'])
end