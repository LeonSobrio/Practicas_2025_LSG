% Task analysis of event-related EEG signals
% Authors: Leon Sobrio, Guiomar Niso
% Cajal Institute (CSIC), Madrid, Spain
% May, 2026

% TOOLBOXES REQUIRED:
%   Signal Processing Toolbox
%   Image Processing Toolbox
%   Computer Vision Toolbox
%   Statistics and ML Toolbox

% This script runs AFTER QualityControl.m and Preprocessing.m on the same protocol.
% It performs event-locked analysis:
%   - Renames raw BIDS triggers to semantic names (user-configurable)
%   - (Optional) Combines stim/response events into compound events
%   - Epochs the data around selected events (user-configurable window)
%   - Bad trials already marked by Preprocessing.m are excluded automatically
%   - Averages by event and/or by user-defined groups
%   - Computes head model (per subject + task)
%   - Computes noise covariance from the baseline of all good epochs
%   - Computes inverse kernel (minimum norm 2018) and applies it to every average
%   - Saves one HTML + one JSON report per (subject, task)

clc; clear;

disp("=== TaskAnalysis script started.")


%% ========================================================================
%                         USER CONFIGURATION PANEL
% =========================================================================

%% ---------------------- PATHS -------------------------------------------

BrainstormDbDir = 'C:\Users\leons\Documents\brainstorm_db';
ReportsDir      = 'C:\Users\leons\Documents\Prácticas\TestingDatasetERP\Task_Reports';
ProtocolName    = 'TestingDatasetERP';   % same protocol as QC and Preprocessing

%% ---------------------- SUBJECTS AND TASKS ------------------------------

% List of participants to analyze (empty = all subjects in protocol)
Subs = {'sub-001', 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-006', ...
        'sub-007', 'sub-008', 'sub-009', 'sub-010', 'sub-011', ...
        'sub-012', 'sub-013', 'sub-014', 'sub-015', 'sub-016', ...
        'sub-017', 'sub-018', 'sub-019', 'sub-020', 'sub-021', ...
        'sub-022', 'sub-023', 'sub-024', 'sub-025', 'sub-026', ...
        'sub-027', 'sub-028', 'sub-029', 'sub-030', 'sub-031', ...
        'sub-032', 'sub-033', 'sub-034', 'sub-035', 'sub-036', ...
        'sub-037', 'sub-038', 'sub-039', 'sub-040'};

% List of BIDS task names to process. The script looks for the processed
% recordings of each task inside each subject's database folder.
TaskNames = {'P3'};

% String that identifies the FINAL processed file in the database (output of
% Preprocessing.m). The script matches it inside the Condition string to
% avoid picking up intermediate files.
ProcessedTag = 'Average_reference';

%% ---------------------- EVENT RENAMING ----------------------------------

% Map raw BIDS trigger names to semantic names. The script renames every
% event listed here before doing anything else.
% Format: Nx2 cell array. Column 1 = original name, column 2 = new name.
% Leave as {} to skip renaming.
%
% Default mapping below corresponds to a go/no-go reversal paradigm.
% Modify or replace entirely for other paradigms.
EventRename = {};

%% ---------------------- COMBINE STIM/RESPONSE ---------------------------

% (Optional) Create new compound events from the succession of two events
% using process_evt_combine. Useful for paradigms where you want to split
% trials by stimulus-response combinations (e.g. correct vs incorrect go).
%
% Format: Nx4 cell array. Each row defines one combination:
%   { name_for_A,  name_for_B,  event_A,  event_B }
% Special keywords for the first two fields:
%   'ignore' = do not create a new event for that side
%   'extend' = create one extended event from A to B
%
% Example (commented):
%   EventCombine = {
%       'ignore', 'go_correct',   'imgC_go_100', 'nCorrect';
%       'ignore', 'go_incorrect', 'imgC_go_100', 'nIncorrect';
%   };
%
% Leave as {} to skip this step.
EventCombine = {};

% Maximum delay between paired events (seconds). Only used if EventCombine
% is not empty.
CombineMaxDelay = 1.0;

%% ---------------------- EPOCHING ----------------------------------------

% Events to epoch around. Comma-separated string (Brainstorm syntax).
% These names must match the renamed events (or the original ones if
% EventRename is empty).
EpochEvents = ['stimulus-11, stimulus-22, stimulus-33, stimulus-44, stimulus-55, ', ...
               'stimulus-12, stimulus-13, stimulus-14, stimulus-15, ', ...
               'stimulus-21, stimulus-23, stimulus-24, stimulus-25, ', ...
               'stimulus-31, stimulus-32, stimulus-34, stimulus-35, ', ...
               'stimulus-41, stimulus-42, stimulus-43, stimulus-45, ', ...
               'stimulus-51, stimulus-52, stimulus-53, stimulus-54'];

% Epoch time window relative to event (seconds)
EpochTime    = [-0.2, 0.5];

% Baseline window for DC offset removal (seconds, relative to event)
BaselineWin  = [-0.2, 0];

%% ---------------------- AVERAGING ---------------------------------------

% Averaging mode:
%   'individual' = one average per event type (one ERP per epoch event)
%   'grouped'    = only averages defined in AvgGroups below
AvgMode = 'both';  %     = both individual averages AND grouped averages
%    AvgMode = 'individual';

% Custom groups of events to average together. Each field is one group;
% the field name becomes the comment of the resulting average.
% Leave as struct() to skip grouped averaging.
%
% Example (commented):
%   AvgGroups.go_all   = {'imgC_go_100', 'imgB_go_90', 'imgE_go_80'};
%   AvgGroups.nogo_all = {'imgfbE_nogo_100', 'imgA_nogo_90', 'imgF_nogo_80'};
%   AvgGroups.new_all  = {'imgfbE_newgo_100', 'imgC_newnogo_100', ...
%                         'imgF_newgo_80', 'imgE_newnogo_80'};

AvgGroups.target    = {'stimulus-11', 'stimulus-22', 'stimulus-33', 'stimulus-44', 'stimulus-55'};
AvgGroups.nontarget = {'stimulus-12', 'stimulus-13', 'stimulus-14', 'stimulus-15', ...
                       'stimulus-21', 'stimulus-23', 'stimulus-24', 'stimulus-25', ...
                       'stimulus-31', 'stimulus-32', 'stimulus-34', 'stimulus-35', ...
                       'stimulus-41', 'stimulus-42', 'stimulus-43', 'stimulus-45', ...
                       'stimulus-51', 'stimulus-52', 'stimulus-53', 'stimulus-54'};

% Average function (for process_average):
%   1 = arithmetic mean   | 2 = absolute value of mean
%   3 = root mean square  | 4 = standard deviation
%   5 = standard error    | 6 = arithmetic mean + standard error
%   7 = median
avg_func = 1;

% Weighted average across trials (1 = yes, 0 = no)
avg_weighted = 0;

%% ---------------------- HEAD MODEL --------------------------------------

% Source space: 1 = cortex surface | 2 = MRI volume | 3 = custom
headModel_sourcespace = 1;

% EEG head model:
%   1 = OpenMEEG BEM  | 2 = 3-shell sphere (Berg)  | 3 = single sphere
headModel_eeg = 2;

% Other modalities (required by process_headmodel; ignored for EEG-only)
headModel_meg   = 3;   % 3 = overlapping spheres
headModel_ecog  = 2;
headModel_seeg  = 2;
headModel_nirs  = 1;

%% ---------------------- NOISE COVARIANCE --------------------------------

noiseCov_sensortypes = 'EEG';
noiseCov_target      = 1;     % 1 = noise covariance | 2 = data covariance
noiseCov_dcoffset    = 1;     % 1 = block-by-block (recommended)
noiseCov_identity    = 0;     % 1 = use identity matrix instead of estimating
noiseCov_copycond    = 0;
noiseCov_copysubj    = 0;
noiseCov_copymatch   = 0;
noiseCov_replacefile = 1;

%% ---------------------- INVERSE SOLUTION (sources) ----------------------

% Inverse method:
%   'minnorm' = minimum norm | 'gls' = generalized least squares | 'lcmv' = beamformer
sources_inverseMethod  = 'minnorm';

% Inverse measure:
%   'amplitude' | 'dspm2018' | 'sloreta' | 'performance'
sources_inverseMeasure = 'amplitude';

% Source orientation:
%   'fixed' = constrained (normal to cortex) | 'free' = unconstrained | 'loose'
sources_sourceOrient   = 'fixed';

sources_loose          = 0.2;
sources_useDepth       = 1;
sources_weightExp      = 0.5;
sources_weightLimit    = 10;
sources_noiseMethod    = 'reg';     % 'reg' | 'diag' | 'none' | 'shrink'
sources_noiseReg       = 0.1;
sources_snrMethod      = 'fixed';   % 'fixed' | 'rms'
sources_snrRms         = 1e-06;
sources_snrFixed       = 3;
sources_computeKernel  = 1;          % 1 = shared kernel (recommended)
sources_dataTypes      = 'EEG';  
sources_snapshotTime   = 0.350;   % Time (s) for source snapshots in the report.
                                  % To view every ms, open the source file in
                                  % the GUI and navigate with the slidebar.

%% ---------------------- FIGURE SIZE -------------------------------------

fig_width  = 1000;
fig_height = 600;
fig_size   = [200, 200, fig_width, fig_height];


%% ========================================================================
%                    DO NOT EDIT BELOW THIS LINE
% =========================================================================

%% Start Brainstorm and load protocol from database
disp('== Start Brainstorm (nogui)');

if ~brainstorm('status')
    brainstorm nogui
end
bst_set('BrainstormDbDir', BrainstormDbDir);
bst_colormaps('RestoreDefaults', 'eeg');

disp(['- BrainstormDbDir:',   bst_get('BrainstormDbDir')]);
disp(['- BrainstormUserDir:', bst_get('BrainstormUserDir')]);
disp(['- HOME env:',          getenv('HOME')]);
disp(['- HOME java:',         char(java.lang.System.getProperty('user.home'))]);

disp(['== Loading protocol: ', ProtocolName]);
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Protocol "', ProtocolName, '" not found in database at: ', BrainstormDbDir, ...
           '. Run QualityControl.m and Preprocessing.m first.']);
end
gui_brainstorm('SetCurrentProtocol', iProtocol);

% Create reports directory if needed
if ~exist(ReportsDir, 'dir'), mkdir(ReportsDir); end


%% Get all data files in protocol
disp('=== Get all data files in protocol');

sFilesAll = bst_process('CallProcess', 'process_select_files_data', [], [], ...
    'subjectname',   '', ...
    'condition',     '', ...
    'tag',           '', ...
    'includebad',    1, ...
    'includeintra',  1, ...
    'includecommon', 0);

if isempty(sFilesAll)
    error('No data files found in protocol. Run QualityControl.m and Preprocessing.m first.');
else
    disp(['Found ', num2str(length(sFilesAll)), ' data files in protocol.']);
end

% Keep only the fully preprocessed files (those tagged with ProcessedTag in
% their condition). These are the output of Preprocessing.m.
condList = {sFilesAll.Condition};
isProcessed = contains(condList, ProcessedTag);
sFilesProc  = sFilesAll(isProcessed);

if isempty(sFilesProc)
    error('No preprocessed files found (tag "%s" not present in any Condition). Run Preprocessing.m first.', ProcessedTag);
end

% Unique participants found in processed files
allParticipants = unique({sFilesProc.SubjectName});


%% ========================================================================
%% Loop over participants
%% ========================================================================
for iSub = 1:length(allParticipants)

    participant = allParticipants{iSub};

    % Skip if not in user-specified subject list
    if ~isempty(Subs) && ~ismember(participant, Subs)
        disp(['=== Skipping participant: ', participant, ' (not in subject list)']);
        continue
    end
    disp(['=== Processing participant: ', participant]);


    %% ====================================================================
    %% Loop over tasks
    %% ====================================================================
    for iTask = 1:length(TaskNames)

        currentTask = TaskNames{iTask};
        disp(['--- Task: ', currentTask, ' (', num2str(iTask), '/', num2str(length(TaskNames)), ')']);

        % Locate the processed file for this (subject, task)
        % We match: subject name + 'task-<currentTask>' + ProcessedTag in Condition
        taskPattern = ['task-', currentTask];
        isThisSubj  = strcmp({sFilesProc.SubjectName}, participant);
        isThisTask  = contains({sFilesProc.Condition}, taskPattern);
        sFilesTask  = sFilesProc(isThisSubj & isThisTask);

        if isempty(sFilesTask)
            warning('No processed file found for %s - %s. Skipping task.', participant, currentTask);
            continue
        end
        if length(sFilesTask) > 1
            warning('Multiple processed files found for %s - %s. Using the first.', participant, currentTask);
            sFilesTask = sFilesTask(1);
        end

        conditionName = sFilesTask.Condition;
        disp(['  Found processed file. Condition: ', conditionName]);

        % Start a fresh report for this (subject, task)
        bst_report('Start', sFilesTask);

        % Initialise JSON structure for this (subject, task)
        jsonData             = struct();
        jsonData.participant = participant;
        jsonData.task        = currentTask;
        jsonData.condition   = conditionName;
        jsonData.protocol    = ProtocolName;
        jsonData.date        = datestr(now, 'yyyy-mm-dd HH:MM:SS');


        % -----------------------------------------------------------------
        %  Rename events
        % -----------------------------------------------------------------
        if ~isempty(EventRename)
            disp('=== Rename events');
            renamedCount = 0;
            for iEvt = 1:size(EventRename, 1)
                oldName = EventRename{iEvt, 1};
                newName = EventRename{iEvt, 2};
                try
                    bst_process('CallProcess', 'process_evt_rename', sFilesTask, [], ...
                        'src',  oldName, ...
                        'dest', newName);
                    renamedCount = renamedCount + 1;
                catch ME
                    warning('Could not rename "%s" -> "%s" for %s: %s', ...
                            oldName, newName, participant, ME.message);
                end
            end
            jsonData.eventsRenamed = renamedCount;
            disp(['  Events renamed: ', num2str(renamedCount), '/', num2str(size(EventRename,1))]);
        else
            jsonData.eventsRenamed = 0;
        end


        % -----------------------------------------------------------------
        %  Combine stim/response (optional)
        % -----------------------------------------------------------------
        if ~isempty(EventCombine)
            disp('=== Combine stim/response events');
            % Build the multi-line string expected by process_evt_combine
            % Each row of EventCombine becomes one line: "f1, f2, f3, f4"
            combineLines = cell(size(EventCombine, 1), 1);
            for iC = 1:size(EventCombine, 1)
                combineLines{iC} = sprintf('%s, %s, %s, %s', ...
                    EventCombine{iC,1}, EventCombine{iC,2}, ...
                    EventCombine{iC,3}, EventCombine{iC,4});
            end
            combineStr = strjoin(combineLines, char(10));

            try
                bst_process('CallProcess', 'process_evt_combine', sFilesTask, [], ...
                    'combine', combineStr, ...
                    'dt',      CombineMaxDelay);
                jsonData.eventCombinations = size(EventCombine, 1);
                disp(['  Combinations applied: ', num2str(size(EventCombine,1))]);
            catch ME
                warning('Combine stim/response failed for %s - %s: %s', ...
                        participant, currentTask, ME.message);
                jsonData.eventCombinations = 0;
            end
        else
            jsonData.eventCombinations = 0;
        end


        % -----------------------------------------------------------------
        %  Convert extended events to simple (point) events
        %  Raw BIDS triggers are stored as extended events (onset + offset).
        %  Brainstorm uses the event duration as epoch window when epoching
        %  extended events, ignoring epochtime. We convert all events to
        %  simple point events (keeping only the onset) before epoching.
        % -----------------------------------------------------------------
        disp('=== Convert extended events to simple (keep onset)');
        allEventNames = EpochEvents;
        bst_process('CallProcess', 'process_evt_simple', sFilesTask, [], ...
            'eventname', allEventNames, ...
            'method',    'start');
        disp('  Extended events converted to point events.');


        % -----------------------------------------------------------------
        %  Epoching: import data segments around events
        %  We pass `condition = currentTask` so the epochs are stored in a
        %  clean folder named after the task (e.g. "inicialfirstrun"),
        %  separated from the raw link folder.
        %  timewindow must be set explicitly to the full raw duration,
        %  otherwise Brainstorm ignores epochtime on BST-DATA raw links.
        % -----------------------------------------------------------------
        disp('=== Read raw file time range');
        rawData    = in_bst_data(sFilesTask.FileName, 'Time');
        rawTimeWin = [rawData.Time(1), rawData.Time(end)];
        disp(['  Raw time window: [', num2str(rawTimeWin(1)), ', ', num2str(rawTimeWin(2)), '] s']);

        disp('=== Import epochs around events');
        sEpochs = bst_process('CallProcess', 'process_import_data_event', sFilesTask, [], ...
            'subjectname',   participant, ...
            'condition',     currentTask, ...
            'eventname',     EpochEvents, ...
            'timewindow',    rawTimeWin, ...
            'epochtime',     EpochTime, ...
            'createcond',    0, ...
            'ignoreshort',   1, ...
            'usectfcomp',    0, ...
            'usessp',        1, ...
            'freq',          [], ...
            'baseline',      BaselineWin, ...
            'blsensortypes', 'EEG');

        if isempty(sEpochs)
            warning('No epochs created for %s - %s. Skipping task.', participant, currentTask);
            continue
        end

        nEpochs = length(sEpochs);
        disp(['  Total epochs created: ', num2str(nEpochs)]);
        jsonData.totalEpochs = nEpochs;


        % -----------------------------------------------------------------
        %  Count good vs bad epochs
        %  Bad segments were already marked by Preprocessing.m on the
        %  continuous data. process_import_data_event inherits those flags,
        %  so epochs that overlap with bad segments are already tagged bad.
        % -----------------------------------------------------------------
        sEpochsGood = bst_process('CallProcess', 'process_select_files_data', [], [], ...
            'subjectname', participant, ...
            'condition',   currentTask, ...
            'tag',         '', ...
            'includebad',  0);

        nGood = length(sEpochsGood);
        nBad  = nEpochs - nGood;
        disp(['  Good epochs: ', num2str(nGood), ' | Bad epochs (inherited from Preprocessing): ', num2str(nBad)]);
        jsonData.goodEpochs = nGood;
        jsonData.badEpochs  = nBad;

        if isempty(sEpochsGood)
            warning('All epochs are marked bad for %s - %s (check Preprocessing bad segments). Skipping task.', ...
                    participant, currentTask);
            continue
        end


        % -----------------------------------------------------------------
        %  Averaging
        % -----------------------------------------------------------------
        disp('=== Compute averages');
        sAvgAll = [];     % accumulates all averages (individual + grouped)
        avgInfo = struct();

        % Individual averages: one per event type
        if strcmpi(AvgMode, 'individual') || strcmpi(AvgMode, 'both')
            disp('  Mode: individual averages (per event)');

            % Group epochs by event name using the file Comment field.
            % Each epoch's Comment starts with the event name.
            epochComments = {sEpochsGood.Comment};
            eventList = unique(cellfun(@(c) extractEventName(c), epochComments, 'UniformOutput', false));

            for iEv = 1:length(eventList)
                evName = eventList{iEv};
                evMask = strcmp(cellfun(@(c) extractEventName(c), epochComments, 'UniformOutput', false), evName);
                sEvent = sEpochsGood(evMask);

                if isempty(sEvent), continue; end

                sAvgEv = bst_process('CallProcess', 'process_average', sEvent, [], ...
                    'avgtype',    1, ...        % 1 = everything together
                    'avg_func',   avg_func, ...
                    'weighted',   avg_weighted, ...
                    'keepevents', 0);

                if ~isempty(sAvgEv)
                    % Tag the average with the event name for clarity
                    bst_process('CallProcess', 'process_add_tag', sAvgEv, [], ...
                        'tag',    ['Avg_', evName], ...
                        'output', 1);  % 1 = add to comment
                    sAvgAll = [sAvgAll, sAvgEv]; %#ok<AGROW>
                    avgInfo.(matlab.lang.makeValidName(['ind_', evName])) = length(sEvent);
                end
            end
        end

        % Grouped averages: combine multiple events into one average
        if (strcmpi(AvgMode, 'grouped') || strcmpi(AvgMode, 'both')) && ~isempty(fieldnames(AvgGroups))
            disp('  Mode: grouped averages');
            groupNames = fieldnames(AvgGroups);
            epochComments = {sEpochsGood.Comment};
            epochEvents   = cellfun(@(c) extractEventName(c), epochComments, 'UniformOutput', false);

            for iG = 1:length(groupNames)
                grpName   = groupNames{iG};
                grpEvents = AvgGroups.(grpName);
                grpMask   = ismember(epochEvents, grpEvents);
                sGroup    = sEpochsGood(grpMask);

                if isempty(sGroup)
                    warning('  No epochs found for group "%s". Skipping.', grpName);
                    continue
                end

                sAvgGrp = bst_process('CallProcess', 'process_average', sGroup, [], ...
                    'avgtype',    1, ...
                    'avg_func',   avg_func, ...
                    'weighted',   avg_weighted, ...
                    'keepevents', 0);

                if ~isempty(sAvgGrp)
                    bst_process('CallProcess', 'process_add_tag', sAvgGrp, [], ...
                        'tag',    ['Group_', grpName], ...
                        'output', 1);
                    sAvgAll = [sAvgAll, sAvgGrp]; %#ok<AGROW>
                    avgInfo.(matlab.lang.makeValidName(['grp_', grpName])) = length(sGroup);
                end
            end
        end

        jsonData.averages = avgInfo;

        if isempty(sAvgAll)
            warning('No averages produced for %s - %s. Skipping source analysis.', participant, currentTask);
            continue
        end

        % Snapshot the first average as a sanity check
        try
            hFigAvg = view_timeseries(sAvgAll(1).FileName, 'EEG');
            set(hFigAvg, 'Position', fig_size);
            bst_report('Snapshot', hFigAvg, sAvgAll(1).FileName, ...
                       ['ERP average example - ', currentTask], fig_size);
            close(hFigAvg);
        catch
            warning('Could not snapshot ERP average for %s - %s.', participant, currentTask);
        end


        % -----------------------------------------------------------------
        %  Head model
        % -----------------------------------------------------------------
        disp('=== Compute head model');
        bst_process('CallProcess', 'process_headmodel', sEpochsGood, [], ...
            'Comment',     '', ...
            'sourcespace', headModel_sourcespace, ...
            'meg',         headModel_meg, ...
            'eeg',         headModel_eeg, ...
            'ecog',        headModel_ecog, ...
            'seeg',        headModel_seeg, ...
            'nirs',        headModel_nirs, ...
            'openmeeg',    struct( ...
                 'BemSelect',    [1, 1, 1], ...
                 'BemCond',      [1, 0.0125, 1], ...
                 'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
                 'BemFiles',     {{}}, ...
                 'isAdjoint',    0, ...
                 'isAdaptative', 1, ...
                 'isSplit',      0, ...
                 'SplitLength',  4000), ...
            'nirstorm',    struct( ...
                 'FluenceFolder',    'https://neuroimage.usc.edu/resources/nst_data/fluence/', ...
                 'smoothing_method', 'geodesic_dist', ...
                 'smoothing_fwhm',   10), ...
            'channelfile', '');


        % -----------------------------------------------------------------
        %  Noise covariance (from baseline of all good epochs)
        % -----------------------------------------------------------------
        disp('=== Compute noise covariance');
        bst_process('CallProcess', 'process_noisecov', sEpochsGood, [], ...
            'baseline',       BaselineWin, ...
            'datatimewindow', EpochTime, ...
            'sensortypes',    noiseCov_sensortypes, ...
            'target',         noiseCov_target, ...
            'dcoffset',       noiseCov_dcoffset, ...
            'identity',       noiseCov_identity, ...
            'copycond',       noiseCov_copycond, ...
            'copysubj',       noiseCov_copysubj, ...
            'copymatch',      noiseCov_copymatch, ...
            'replacefile',    noiseCov_replacefile);

        % Visualize and snapshot the noise covariance
        try
            [sStudy, ~] = bst_get('AnyFile', sEpochsGood(1).FileName);
            if ~isempty(sStudy) && ~isempty(sStudy.NoiseCov)
                NoiseCovFile = sStudy.NoiseCov(1).FileName;
                hFigNoise    = view_noisecov(NoiseCovFile);
                set(hFigNoise, 'Position', fig_size);
                bst_report('Snapshot', hFigNoise, NoiseCovFile, ...
                           ['Noise covariance - ', currentTask], fig_size);
                close(hFigNoise);
            end
        catch
            warning('Could not snapshot noise covariance for %s - %s.', participant, currentTask);
        end


        % -----------------------------------------------------------------
        %  Inverse kernel (shared across all averages of this task)
        %  IMPORTANT: process_inverse_2018 expects the averaged files
        %  (sAvgAll), not the raw epochs. With 'output = 1' (shared kernel),
        %  it returns the source files corresponding to each average.
        % -----------------------------------------------------------------
        disp('=== Compute inverse kernel (sources)');
        sSources = bst_process('CallProcess', 'process_inverse_2018', sAvgAll, [], ...
            'output',  1, ...  % 1 = kernel only (shared across files)
            'inverse', struct( ...
                 'Comment',        'MN: EEG', ...
                 'InverseMethod',  sources_inverseMethod, ...
                 'InverseMeasure', sources_inverseMeasure, ...
                 'SourceOrient',   {{sources_sourceOrient}}, ...
                 'Loose',          sources_loose, ...
                 'UseDepth',       sources_useDepth, ...
                 'WeightExp',      sources_weightExp, ...
                 'WeightLimit',    sources_weightLimit, ...
                 'NoiseMethod',    sources_noiseMethod, ...
                 'NoiseReg',       sources_noiseReg, ...
                 'SnrMethod',      sources_snrMethod, ...
                 'SnrRms',         sources_snrRms, ...
                 'SnrFixed',       sources_snrFixed, ...
                 'ComputeKernel',  sources_computeKernel, ...
                 'DataTypes',      {{sources_dataTypes}}));

        if isempty(sSources)
            warning('Inverse computation failed for %s - %s.', participant, currentTask);
            continue
        end


        % -----------------------------------------------------------------
        %  Visualize sources on the first average (sanity check)
        % -----------------------------------------------------------------
        disp('=== Visualize sources (first average)');
        try
            % sSources(1).FileName is already a valid Brainstorm source link
            hFigSrc = view_surface_data([], sSources(1).FileName);
            panel_time('SetCurrentTime', sources_snapshotTime);
            figure_3d('SetStandardView', hFigSrc, 'left');
            bst_report('Snapshot', hFigSrc, sSources(1).FileName, ...
                       ['Sources Left - ', currentTask], fig_size);
            figure_3d('SetStandardView', hFigSrc, 'right');
            bst_report('Snapshot', hFigSrc, sSources(1).FileName, ...
                       ['Sources Right - ', currentTask], fig_size);
            figure_3d('SetStandardView', hFigSrc, 'top');
            bst_report('Snapshot', hFigSrc, sSources(1).FileName, ...
                       ['Sources Top - ', currentTask], fig_size);
            close(hFigSrc);
        catch ME
            warning('Source visualization failed for %s - %s: %s', ...
                    participant, currentTask, ME.message);
        end


        % -----------------------------------------------------------------
        %  Save HTML + JSON reports for this (subject, task)
        % -----------------------------------------------------------------
        disp('=== Save HTML report');
        htmlName   = fullfile(ReportsDir, sprintf('Task-%s-%s-%s.html', ...
                              participant, conditionName, ProtocolName));
        ReportFile = bst_report('Save', []);
        bst_report('Export', ReportFile, htmlName);
        disp(['  HTML saved: ', htmlName]);

        disp('=== Save JSON report');
        jsonFile = fullfile(ReportsDir, sprintf('Task-%s-%s-%s.json', ...
                            participant, conditionName, ProtocolName));
        fid = fopen(jsonFile, 'w');
        if fid ~= -1
            fprintf(fid, '%s', jsonencode(jsonData, PrettyPrint=true));
            fclose(fid);
            disp(['  JSON saved: ', jsonFile]);
        else
            warning('Could not open JSON file for writing: %s', jsonFile);
        end

    end  % iTask

    disp(['=== Finished participant: ', participant]);

end  % iSub

disp('=== TaskAnalysis complete. Pipeline finished.');


%% ========================================================================
%                         HELPER FUNCTIONS
% =========================================================================

function evName = extractEventName(comment)
    % Extract the event name from a Brainstorm epoch Comment string.
    % Epoch comments are typically of the form:
    %   "<eventname> (#NN)" or "<eventname>"
    % We strip trailing parentheses and surrounding whitespace.
    evName = regexprep(comment, '\s*\(.*\)\s*$', '');
    evName = strtrim(evName);
end