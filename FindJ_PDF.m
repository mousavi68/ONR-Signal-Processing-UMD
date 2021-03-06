function [J_Phase2,flag] = FindJ_PDF(J_PEP,J_A,xi,f)
% Description:----------------------------------------
% Function to find BCG/SCG_HF J-peak that has highest PDF

% Inputs: -----------------------------------------------
% J_PEP: timing b/w ECG R-peak and J-wave
% J_A: amplitude of J-waves
% f: probability density estimate
% xi: vector of sample data

% Outputs:---------------------------------------------
% J_Phase2: final acceptable J-wave selected after implementing PDF based method
% flag: has three values: -1,0,1 --> 1 means current beat has an acceptable
% J-wave with PDF>0, 0 means none of the J-wave candidates had PDF>0, -1
% means from the start no J-peak candidates were extracted for this beat 

% Body of the code:----------------------------------
for ii = 1:length(J_PEP)
    if sum(~isnan(J_PEP(:,ii)))>0
        [temp,tempIdx] = min(abs(transpose(J_PEP(:,ii)-xi)));   
        dd = tempIdx(~isnan(temp)); % Drop indices with Nan
        [pdf,idx] = max(f(dd)); % find the J-wave with highest PDF
        J_Phase2.Val(ii,:) = [J_PEP(idx,ii),J_A(idx,ii),idx,pdf];
        if pdf >0
            flag(ii,1) = 1; % One acceptable J-wave with PDF>0 is selected
        else
            J_Phase2.Val(ii,:) = nan(1,4);
            flag(ii,1) = 0; % all J-wave candidates in this beats have zero PDF
        end
    else
        J_Phase2.Val(ii,:) = nan(1,4);
        flag(ii,1) = -1;   % didn't detect any J-peaks in the beat from the start
    end
    J_Phase2.Label = {'J_PEP','J_amp','Index','PDF'};
end
end