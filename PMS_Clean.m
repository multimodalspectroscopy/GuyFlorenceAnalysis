%% PMS_Clean - clean PMS datasets by replacing anomalous values with NaN
% This script imports the PMS Clinical Datasets containing synchronised 
% timetables with anomalous NIRS, DCS, and systemic data
% The anomalous data exhibits a periodic pattern of max absorption in all
% 11 variables within each file
% This script removes these periodic anomalous areas and replaces them with
% NaN values to clean the data

% Cleansed PMS files are saved in both .mat and .xlsx formats, following
% same nomenclature as previous.

% N.B. plots used in for-loop. This was done to aid script development and
% to visualise methods and distributions. Remove plots before executing to
% speed up.
% Additional TEST plots can be used to ensure data has been edited properly
    % See plotFLORENCE_PMS function.

% Author: Guy Murphy
% Date of latest version: 27/03/24

%% Define Folder to Scan and Save Location

% Define folder to scan containing uncleaned PMS data
FolderToScan = '/Users/guymurphy/Documents/MATLAB/MPHY0012_ResearchProject/Data/Clinical Datasets/PMS'; % PMS data folder

% Define folder to save the cleaned data - 'PMS_Cleaned'
SaveDirectoryCleanedData = '/Users/guymurphy/Documents/MATLAB/MPHY0012_ResearchProject/Data/Clinical Datasets/PMS/PMS_Cleaned'; % save folder
if ~exist(SaveDirectoryCleanedData,'dir')
    mkdir(SaveDirectoryCleanedData)
end

cd(FolderToScan)
PMSfileList = dir(fullfile(FolderToScan,'*.mat'));
NbSubject = length(PMSfileList);

%% ****START: Import & Cleaning****

for i = 21:NbSubject %1:NbSubject or indicate which file you want to process (NbSubject)
    % N.B. use single file indication when testing code and full NbSubject
    % from PMSfileList when running

    %% Import data
    data = load(fullfile(FolderToScan, PMSfileList(i).name));
        % Each 'PMS_0**B_Session*_ClinicalDataset' contains a single
        % variable; a timetable containing the 11 data variables (12 for
        % PMS_001B).
        % For each of the 23 subjects, there is only 1 session.

    % Extract dataset name + data to access dataset in struct
    datasetName = fieldnames(data);
    dataset = data.(datasetName{1});

    % Extract subject + session number for plots and filenames
    expression = 'PMS_(\d+B)_Session(\d+)';
    matches = regexp(datasetName, expression,'tokens');
    subjectNum = matches{1}{1}{1};
    sessionNum = matches{1}{1}{2};
    
    dataInfo.subject = subjectNum;
    dataInfo.session = sessionNum;

    %% TEST PLOT - before cleaning
    plotFLORENCE_PMS(dataset,dataInfo.subject,dataInfo.session,1,2); %remove before executing

    %% Identify and replace anomalous values
    % Anomalous values are at same time points across all metrics

    % Identify anomalous values in [oxCCO]_Changes by replacing all values
    % where the change in [oxCCO]_Changes from the preceding point is 
    % greater than a certain threshold (diffThreshold), as well also all 
    % values with magnitude greater than another certain threshold 
    % (magThreshold).

    % Index and replace all values at these time points with 'NaN' in other
    % metrics.

    % [oxCCO]_Changes most suitable metric to base exclusion criteria on as
    % it typically remains fairly stable with smaller changes in the good
    % data regions and much larger changes in the anomalous data ranges.
    % Therefore, we can use lower automatic detection threshold values.

    % Define threshold values
    % Default values: diffThreshold = 5; magThreshold = 20;
    diffThreshold = 5;
    magThreshold = 20;

    % Manually altered threshold values for specific subjects/sessions
    if i == 2
        diffThreshold = 2.5;
    elseif i == 3
        diffThreshold = 3.5;
    elseif i == 6
        diffThreshold = 2.5;
        magThreshold = 40;
    elseif i == 7
        diffThreshold = 2.5;
    elseif i == 8
        diffThreshold = 60;
        magThreshold = -200;
    elseif i == 9
        diffThreshold = 40;
        magThreshold = 400;
    elseif i == 10
        diffThreshold = 3.5;
    elseif i == 11
        diffThreshold = 7.5;
        magThreshold = 60;
    elseif i == 16
        diffThreshold = 2.5;
    elseif i == 17
        diffThreshold = 3;
    elseif i == 18
        diffThreshold = 15;
        magThreshold = -20;
    elseif i == 19
        diffThreshold = 7.5;
    elseif i == 20
        diffThreshold = 15;
        magThreshold = 100;
    elseif i == 21
        diffThreshold = 10;
        magThreshold = 47.5;
    elseif i == 22
        diffThreshold = 1e13;
        magThreshold = 1.5e14;
    elseif i == 23
        magThreshold = 30;
    end
    
    % Clean datasets
    if i == 8 || i == 18
        condition = abs(diff(dataset.("[oxCCO]_Changes"))) > diffThreshold; %define condition
        dataset{condition, 1:end} = NaN;

        % Remove absolute condition for dataset 8 and 18 (PMS_018B +
        % PMS_043B - negative maximum values)
        condition = dataset.("[oxCCO]_Changes") > magThreshold; %define condition
        dataset{condition, 1:end} = NaN;

    else
        % Default conditions
        % Replace dataset values with NaN at time points where the absolute
        % magnitude of the difference between two sequential values in
        % [oxCCO]_Changes > |diffThreshold|, for all variables in dataset
        condition = abs(diff(dataset.("[oxCCO]_Changes"))) > diffThreshold; %define condition
        dataset{condition, 1:end} = NaN;

        % Replace dataset values with NaN at time points where the 
        % absolute magnitude of [oxCCO]_Changes > |magThreshold|, for all
        % variables in dataset
        condition = abs(dataset.("[oxCCO]_Changes")) > magThreshold; %define condition
        dataset{condition, 1:end} = NaN;

    end

    % ****END: Import & Cleaning****

    %% TEST PLOT - after cleaning
    plotFLORENCE_PMS(dataset,dataInfo.subject,dataInfo.session,3,4); %remove before executing

    %% Save ClinicalDataset_CLEAN
    
    cd(SaveDirectoryCleanedData)

    % SAVE TO EXCEL
    excelFileName = ['PMS_' dataInfo.subject '_Session' dataInfo.session '_ClinicalDatasetCleaned.xlsx'];
    writetimetable(dataset,fullfile(SaveDirectoryCleanedData,excelFileName))

    % SAVE TO .MAT
        % Assign Patient Name + Session to cleaned dataset timetable
    matVarName = strcat('PMS_',dataInfo.subject,'_Session',dataInfo.session);
    eval([matVarName ' = dataset;']);

        % Assign File Name
    matFileName = ['PMS_' dataInfo.subject '_Session' dataInfo.session '_ClinicalDatasetCleaned.mat'];
    
        % Save cleaned dataset variable as a .mat file
    save(fullfile(SaveDirectoryCleanedData,matFileName), matVarName)
    
end