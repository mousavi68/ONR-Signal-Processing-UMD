function param = SetParameters(fs,SID)
% Description:-----------------------------------------
% This function evaluates parameters needed in the main code for feature
% extraction. Parameters evaluated in this function include: threshold values for motion artifact, property values for the findpeaks function

% Inputs:------------------------------------------------
% fs: sampling frequency
% SID: current subject identification #

% Outputs:---------------------------------------------
% param: structure containing the variables needed in the analysis of
% motion artifact and feature extraction, the parameters are defined for
% each subject separately using if/elesif coding

% Body of the code:----------------------------------
%%% Subject 010-----------------------------------
if SID == 10
    % Motion Artifact
    param.MotMult = 10; 
    param.MotMultSCG = 10; 
    % ECG
    param.ECG.RMinPeakHeight = 0.02;
    % BCG
    param.BCG_A.JMinPeakHeight = 5e-5;
    param.BCG_A.IMinPeakHeight = 1e-4;
    param.BCG_A.KMinPeakHeight = 1e-4;
    param.BCG_A.HMinPeakHeight = 1e-4;
    param.BCG_A.LMinPeakHeight = 1e-4;
    param.BCG_A.JSt = 0.05*fs;
    param.BCG_A.ISt = 0.02*fs;
    param.BCG_A.HSt = 0.01*fs;
    param.BCG_A.JEnd = 0.2*fs;
    param.BCG_A.KEnd = 0.4*fs;
    param.BCG_A.LEnd = 0.1*fs;
    param.BCG_A.JProm = 1e-4;
    % SCG
    param.SCG.JMinPeakHeight = 1e-4;
    param.SCG.IMinPeakHeight = 1e-4;
    param.SCG.KMinPeakHeight = 1e-4;
    param.SCG.HMinPeakHeight = 1e-4;
    param.SCG.LMinPeakHeight = 1e-4;
    param.SCG.JSt = 0.1*fs;
    param.SCG.ISt = 0.02*fs;
    param.SCG.HSt = 0.01*fs;
    param.SCG.JEnd = 0.4*fs;
    param.SCG.KEnd = 0.4*fs;
    param.SCG.LEnd = 0.1*fs;
    param.SCG.JProm = 1e-4;
    param.AOSt = 0.04*fs;
    param.AOEnd = 0.1*fs;
    param.AOMinPeakHeight = 5e-5;
    param.ACSt = 0.2*fs;
    param.ACEnd = 0.5*fs;
    param.ACMinPeakHeight = 5e-5;
    % BP related parameters
    param.BP.SPMinPeakHeight = 85;
    param.BP.DPMinPeakHeight = -120;
    param.BP.SPEnd = 0.5*fs;
    param.BP.SPSt = 0.1*fs;
    % PPG Related parameters
    param.PPG.PPGMinPeakHeight = 0.15;
    param.PPG.PPGFtMinPeakHeight = 0.05;
    param.PPG.PPGPkEnd = 0.5*fs;
    param.PPG.PPGPkSt = 0.1*fs;
%%% subject 011-----------------------------------
elseif SID == 11
    % Motion Artifact
    param.MotMult = 10; 
    param.MotMultSCG = 10; 
    % ECG
    param.ECG.RMinPeakHeight = 0.02;
    % BCG
    param.BCG_A.JMinPeakHeight = 5e-5;
    param.BCG_A.IMinPeakHeight = 1e-4;
    param.BCG_A.KMinPeakHeight = 1e-4;
    param.BCG_A.HMinPeakHeight = 1e-4;
    param.BCG_A.LMinPeakHeight = 1e-4;
    param.BCG_A.JSt = 0.05*fs;
    param.BCG_A.ISt = 0.02*fs;
    param.BCG_A.HSt = 0.01*fs;
    param.BCG_A.JEnd = 0.2*fs;
    param.BCG_A.KEnd = 0.4*fs;
    param.BCG_A.LEnd = 0.1*fs;
    param.BCG_A.JProm = 1e-4;
    % SCG
    param.SCG.JMinPeakHeight = 1e-4;
    param.SCG.IMinPeakHeight = 1e-4;
    param.SCG.KMinPeakHeight = 1e-4;
    param.SCG.HMinPeakHeight = 1e-4;
    param.SCG.LMinPeakHeight = 1e-4;
    param.SCG.JSt = 0.1*fs;
    param.SCG.ISt = 0.02*fs;
    param.SCG.HSt = 0.01*fs;
    param.SCG.JEnd = 0.4*fs;
    param.SCG.KEnd = 0.4*fs;
    param.SCG.LEnd = 0.1*fs;
    param.SCG.JProm = 1e-4;
    param.AOSt = 0.03*fs;
    param.AOEnd = 0.1*fs;
    param.AOMinPeakHeight = 5e-5;
    param.ACSt = 0.2*fs;
    param.ACEnd = 0.5*fs;
    param.ACMinPeakHeight = 5e-5;
    % BP related parameters
    param.BP.SPMinPeakHeight = 120;
    param.BP.DPMinPeakHeight = -120;
    param.BP.SPEnd = 0.5*fs;
    param.BP.SPSt = 0.1*fs;
    % PPG Related parameters
    param.PPG.PPGMinPeakHeight = 0.15;
    param.PPG.PPGFtMinPeakHeight = 0.05;
    param.PPG.PPGPkEnd = 0.5*fs;
    param.PPG.PPGPkSt = 0.1*fs;
end
end