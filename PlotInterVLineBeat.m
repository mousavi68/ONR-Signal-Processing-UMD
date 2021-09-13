function [] = PlotInterVLineBeat(ylim,timesig,subjectInterventionXLS,fileIndex,eventData,startIdxCut,fs,textloc,textOn,startBeat,endBeat)
% Description:----------------------------------------
% This function is used to put vertical lines on the plot with Beat number as X-axis to mark the start
% of each intervention. The function can also print the name of each
% intervention on the top of the figure

% Inputs:-----------------------------------------------
% ylim: upper limit of the current axis
% timesig: signal in the x-axis
% subjectInterventionXLS: name of the excel file containing intervention timing information 
% fileIndex: current subject index in the excel file
% eventData: variable containing subject information
% startIdxCut: vector containing intervention start sample # 
% fs: sampling frequency
% textOn: flag to not to print text or not to print (set 1 for print, set 0 for not print)
% textloc: controls the vertical location of the text 

% Outputs:---------------------------------------------
% no outputs

%Body of the Code----------------------------------
eventIdx = 1;
validEventIdx = 0;
while ~isnan(subjectInterventionXLS{fileIndex+1,2+(eventIdx-1)*4}) %eventIdx<=4 %
    name = ['event' num2str(eventIdx)];
    if eventData.(name).note == 'Y'
            validEventIdx = validEventIdx+1;
            xline(startBeat(validEventIdx))
            if textOn == 1
                if validEventIdx<length(startIdxCut)
                    text((startBeat(validEventIdx)+startBeat(validEventIdx+1))/2,ylim*textloc,eventData.(name).eventName)
                else
                    text((startBeat(validEventIdx)+endBeat)/2,ylim*textloc,eventData.(name).eventName)
                end
            end
    end
     eventIdx = eventIdx+1;
end