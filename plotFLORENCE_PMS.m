function [] = plotFLORENCE_PMS(dataset,subject,session,figNum1, figNum2)
% plotFLORENCE - provides plots of PMS datasets to visually inspect which
%                 are anomalous, with subject and session number.
%                 N.B. PMS dataset is considered anomalous if a periodic
%                 pattern of max absorption signal is detected.

% INPUT VARIABLES:
    % dataset:  - Synchronised timetable of FLORENCE data
    % subject:  - Subject number for recording (for plot titles)
    % session:  - Session number for recording
    % figNum1/2: - Manual figure input to allow for repeated use of
    %              function in a single script

% PLOTS:
%   - Changes in Chromophore Concentrations (HbO2 + HHb + oxCCO)
%       + Changes in Chromophore Derivations (HbDiff + HbT)
%       + Absolute Chromophore Concentrations (Hb02 + HHb)
%   - Changes in BFi, StO2, Temperature

%% Define Variables

% Extract time and variables from dataset
Time = dataset.Time;
HbO2_delta = dataset.("[HbO2]_Changes");
HHb_delta = dataset.("[HHb]_Changes");
oxCCO_delta = dataset.("[oxCCO]_Changes");
HBT_delta = dataset.("[HBT]_Changes");
HBdiff_delta = dataset.("[HBdiff]_Changes");

HbO2_abs = dataset.("[HbO2]_Absolute");
HHb_abs = dataset.("[HHb]_Absolute");
HBT_abs = dataset.("[HBT]_Absolute");

BFi = dataset.BFi;
StO2 = dataset.StO_2;
Temp = dataset.Temperature;


%% Plots

% Chromophore Concentrations (change + absolute)
figure(figNum1)
sgtitle(['Subject: ' subject ', Session: ' session])
subplot(3,1,1)
plot(Time,HbO2_delta,'r', Time,HHb_delta,'b', Time,oxCCO_delta,'g')
title("Changes in Chromophore Concentrations")
xlabel("Time")
legend('HbO2','HHb','oxCCO','Location','southwest')

subplot(3,1,2)
plot(Time,HBdiff_delta,'b', Time,HBT_delta,'r')
title('Changes in HBdiff and HBT')
xlabel("Time")
legend('HbDiff','HbT','Location','southwest')

subplot(3,1,3)
plot(Time,HbO2_abs,'r', Time,HHb_abs,'b', Time,HBT_abs,'g')
title("Absolute Chromophore Concentrations")
xlabel("Time")
legend('HbO2','HHb','HbT','Location','southwest')

% Bloodflow (DCS), StO2, Temperature
figure(figNum2)
sgtitle(['Subject: ' subject ', Session: ' session])
subplot(3,1,1)
plot(Time,BFi,'m')
title("Bloodflow Index (BFi)")
xlabel("Time")

subplot(3,1,2)
plot(Time,StO2,'r')
title("Absolute Tissue Oxygen Saturation (StO2)")
xlabel("Time")

subplot(3,1,3)
plot(Time,Temp,'c')
title("Temperature")
xlabel("Time")

end