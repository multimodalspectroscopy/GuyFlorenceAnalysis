%% FLORENCE_DistributionAnalysis - Analyse distribution of absolute tissue oxygenation
% This script imports the BBS and PMS_Cleaned Clinical Datasets containing
% synchronised timetables of NIRS, DCS, and systemic data, with the major
% light source artefacts removed for PMS data.
% It then uses various statistical tests to fit the data and test for
% normality. It identifies num of modes. Size and magnitude of multiple
% modes, and degree of separation from normality, correlate with the
% stability of the recording, i.e. uneventful/eventful.

% Author: Guy Murphy
% Date of latest version: 10/04/24

% N.B. plots used in for-loop. This was done to aid script development and
% to visualise methods and distributions. Remove plots before executing to
% speed up.
% Additional TEST plots can be used to ensure data has been edited properly
    % See plotFLORENCE function.

%% Define Folder to Scan

% Define folder to scan containing BBS and PMS_Cleaned datasets
FolderToScan = '/Users/guymurphy/Documents/MATLAB/MPHY0012_ResearchProject/Data/Clinical Datasets/All_BBS_PMSclean';
    % N.B. BBS_036_Session1 removed as Filesize = 0b
ClinDataDir = '/Users/guymurphy/Documents/MATLAB/MPHY0012_ResearchProject/Data/Clinical Datasets';

fileList = dir(fullfile(FolderToScan,'*.mat'));
NbSubject = length(fileList);

%% Pre-allocate arrays for table creation
    % Slower, but easier to manage data using for-loop and pre-allocating

% Pre-allocate array for dataset names (only needed when creating table at end)
recording = cell(NbSubject,1);

% Pre-allocate arrays for goodness-of-fit test results
adDistTest = ones(NbSubject,1); % Anderson-Darling, alpha = 0.05 (default)
adDistTest2 = ones(NbSubject,1); % Anderson-Darling, alpha = 0.0005

lillieDistTest = ones(NbSubject,1); % Lilliefors, alpha = 0.05 (default)
lillieDistTest2 = ones(NbSubject,1); % Lilliefors, alpha = 0.0005

jbDistTest = ones(NbSubject,1); % Jarque-Bera, alpha = 0.05 (default)
jbDistTest2 = ones(NbSubject,1); % Jarque-Bera, alpha = 0.0005

% Pre-allocate array for num_modes result
numModesList = zeros(NbSubject,1);


%% ****START****

for i = 1:NbSubject %1:NbSubject or indicate which file you want to process (NbSubject)
    % N.B. use single file indication when testing code and full NbSubject when running

    %% Import data
    data = load(fullfile(FolderToScan, fileList(i).name));

    % Extract dataset name + data to access dataset in struct
    datasetName = fieldnames(data);
    dataset = data.(datasetName{1});
    
    % Create list of dataset names for table
    recording(i) = datasetName;

    %% Extract study name and subject + session number for plots and filenames

    BBSexpression = '(\w\w\w)_(\d+)_Session(\d+)';
    PMSexpression = '(\w\w\w)_(\d+B)_Session(\d+)';
    BBSmatches = regexp(datasetName, BBSexpression,'tokens');
    PMSmatches = regexp(datasetName, PMSexpression,'tokens');

    if ~isempty(BBSmatches{1})
        dataInfo.study = BBSmatches{1}{1}{1};
        dataInfo.subject = BBSmatches{1}{1}{2};
        dataInfo.session = BBSmatches{1}{1}{3};
    elseif ~isempty(PMSmatches{1})
        dataInfo.study = PMSmatches{1}{1}{1};
        dataInfo.subject = PMSmatches{1}{1}{2};
        dataInfo.session = PMSmatches{1}{1}{3};
    end
    
    %% TEST PLOT - before removal of out-of-range values
    plotFLORENCE(dataset,dataInfo.study,dataInfo.subject,dataInfo.session,1,2);

    %% Remove out of range values

    % Logical indexing to remove values outside expected range for StO_2
        % Expected range: 40-90 
    StO2_lowthreshold = 40;
    StO2_highthreshold = 90;

    StO2_index = (dataset.StO_2 >= StO2_lowthreshold & dataset.StO_2 <= StO2_highthreshold);
        % 1/true for in range, 0/false for out of range

    dataset.StO2_index = StO2_index;

    % Save indexed, in range (IR) StO2 data in new variable
    StO2_IR = dataset.StO_2(dataset.StO2_index==1);

    %% TEST PLOT - after removal of out-of-range values
    figure(3);
    clf(3);
    plot(dataset.Time(dataset.StO2_index==1),StO2_IR)

    %% (1) Visual Inspection: Statistical Analysis for Absolute Tissue Oxygen Saturation (StO_2)
    
    if ~isempty(StO2_IR) % some StO2 datasets have no data in range meaning cannot fit normal
                         % e.g. BBS_066_Session1+2
        
        % Normal distribution fit to StO2_IR to estimate mu and sigma
        StO2_IRNormalFit = fitdist(StO2_IR,'Normal');

        % Plot empirical StO2 cdf and hypothesised StO2 cdf using mu and sigma
        % from norm dist fit
        x = linspace(min(StO2_IR),max(StO2_IR),1000);
        y = normcdf(x,StO2_IRNormalFit.mu,StO2_IRNormalFit.sigma); % hypothesised StO2 cdf

        figure(4);
        clf(4); % clear fig each loop iteration to keep 'color' the same for cdfplot
        cdfplot(StO2_IR), hold on; % empirical StO2 cdf
        title(['StO2 CDFs: ' dataInfo.study '\_' dataInfo.subject '\_Session' dataInfo.session])
        plot(x,y,'r'), hold off; % hypothesised StO2 cdf
        legend('Empirical StO2 CDF','Hypothesised StO2 CDF')
    

        % Plot StO2_IR histogram + normal distribution fit
        figure(5);
        clf(5); % clear fig each loop iteration to keep 'color' the same for histogram
        plot(StO2_IRNormalFit)
        title(['StO2 Distribution Fits: ' dataInfo.study '\_' dataInfo.subject '\_Session' dataInfo.session])
        xlabel('StO2 Data Values'), hold on;

        % Plot KDE non-parametric fit to StO2_IR histogram
            % Kernel density estimation
        [StO2_f, StO2_xi] = ksdensity(StO2_IR);
            % Apply kernel density/non-parametric fit over hist + norm dist
        plot(StO2_xi,StO2_f,'Color','k')
        legend('StO2 Histogram','Normal Distribution Fit','Kernel Density Estimation')


        %% (2) Quantitative Inspection: Statistical Analysis for Absolute Tissue Oxygen Saturation (StO_2)
    
        % - Goodness-of-fit hypothesis testing, e.g. Anderson-Darling, Lilliefors, Jarque-Bera Test
        % - Distribution tests for multimodality
        %
        % Assumption: with no artefacts or desaturation events, the StO2 data
        % for each recording should be normally distributed. If there are
        % artefacts in the data or desaturation events, there should be a
        % multimodal distribution of the data.
        %
        % Therefore, testing for normality is a proposed method of quickly
        % identifying eventful/uneventful recording sessions.

        %% (2a) Goodness-Of-Fit Tests for Normal

        % Anderson-Darling (adtest), Lilliefors (lillietest), & Jarque-Bera (jbtest) tests used

        % 0 => Null Hypothesis (H0): StO2 dataset is normally distributed
        % 1 => Alternative Hypothesis (H1): StO2 dataset is not normally
        % distributed

        % Default significance level (alpha): 0.05 (5%)

        % Anderson-Darling
            % default alpha
        [ad_H,ad_P,ad_Stat,ad_CV] = adtest(StO2_IR,'Distribution','norm');
        adDistTest(i) = ad_H;
            % alpha = 0.05%
        [ad2_H,ad2_P,ad2_Stat,ad2_CV] = adtest(StO2_IR,'Distribution','norm','Alpha',0.0005);
        adDistTest2(i) = ad2_H;

        % Lilliefors
            % default alpha
        [lillie_H,lillie_P,lillie_Stat,lillie_CV] = lillietest(StO2_IR,'Distr','norm');
        lillieDistTest(i) = lillie_H;
            % alpha  = 0.1% - lowest valid sig level for lillietest
        [lillie2_H,lillie2_P,lillie2_Stat,lillie2_CV] = lillietest(StO2_IR,'Distr','norm','Alpha',0.001);
        lillieDistTest2(i) = lillie2_H; 

        % Jarque-Bera
            % default alpha
        [jb_H,jb_P,jb_Stat,jb_CV] = jbtest(StO2_IR);
        jbDistTest(i) = jb_H;
            % alpha = 0.1% - lowest valid sig level for jbtest
        [jb2_H,jb2_P,jb2_Stat,jb2_CV] = jbtest(StO2_IR,0.001);
        jbDistTest2(i) = jb2_H;            

        %% (2b) Test for Multimodality
            % Mode counter

        % Identify number of modes (i.e 'peaks') in each dataset using 
        % gradient changes in kernel distribution estimation (see line ***118***)
        dStO2_f = gradient(StO2_f,StO2_xi);
        f_threshold = 0.015; %threshold for gradient changes (to remove small ones from data noise)
        dStO2_f_thresholded = dStO2_f .* (StO2_f >= f_threshold);

        negative_changes = diff(sign(dStO2_f_thresholded)) < 0; % specify only -ve changes in gradient, i.e. peaks
        num_modes = nnz(negative_changes);
        numModesList(i)= num_modes;

    end
end
% ****END: Import & Analysis****

%% Visual Inspection Results - ignore/remove
% Normalilty test by visual inspection
% Results added manually
% 0 (H0) for normally distributed
% 1 (H1) for not normally distributed (reject null hypothesis)

% Pre-allocate
visInspTest = ones(NbSubject,1);

% Asign 0 values at indices of datasets that visually appear normally distributed
visNorm_indices = [1,2,3,4,7,8,11,12,13,27,31,33,39,59,61,62,71,73,78,80,...
    84,86,89,90,91,93,95,101,102,104,105,159,164];
visInspTest(visNorm_indices) = 0;


%% Save results to table

% Create table variable for subject number (i) used for each dataset in code
subject_numbers = (1:NbSubject)';

% Generate table of analysis data
DistributionData = table(subject_numbers,recording,numModesList,adDistTest, ...
    adDistTest2,lillieDistTest,lillieDistTest2,jbDistTest,jbDistTest2,visInspTest);

% Rename table variables
varNames = DistributionData.Properties.VariableNames;
newVarNames = {'Session Number in Script (i)','Session ID','Number of Modes',...
    'AD Test alpha=0.05','AD alpha=0.0005','Lilliefors Test alpha=0.05',...
    'Lilliefors alpha=0.001','JB Test alpha=0.05','JB alpha=0.001','Norm. Dist. Visual ID'};

for j = 1:length(varNames)
    DistributionData = renamevars(DistributionData,varNames{j},newVarNames{j});
end

% SAVE TO EXCEL
cd(ClinDataDir)
excelFileName = 'FLORENCE_DistributionAnalysisResults.xlsx';
writetable(DistributionData,fullfile(ClinDataDir,excelFileName))
