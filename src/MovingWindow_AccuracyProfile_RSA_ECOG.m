function M = MovingWindow_AccuracyProfile_RSA_ECOG(condition, subjects, holdouts, varargin)
% NOTE: If running as a batch of jobs, the set of windows should be held
% constant over the batch, unless you define sub-directories where chuncks
% of windows will be written.
% DEPENDENCIES
% 1. glmnet
% 2. sqrt_truncate_r.m (from WholeBrain_RSA/src)
    p = inputParser();
    if isdeployed
        addRequired(p, 'condition');
        addOptional(p, 'subjects', '1 2 3 5 7 8 9 10', @ischar);
        addOptional(p, 'holdouts', '1 2 3 4 5 6 7 8 9 10', @ischar);
        addParameter(p, 'windows', '10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500 510 520 530 540 550 560 570 580 590 600 610 620 630 640 650 660 670 680 690 700 710 720 730 740 750 760 770 780 790 800 810 820 830 840 850 860 870 880 890 900 910 920 930 940 950 960 970 980 990 1000 1010 1020 1030 1040 1050 1060 1070 1080 1090 1100 1110 1120 1130 1140 1150 1160 1170 1180 1190 1200 1210 1220 1230 1240 1250 1260 1270 1280 1290 1300 1310 1320 1330 1340 1350 1360 1370 1380 1390 1400 1410 1420 1430 1440 1450 1460 1470 1480 1490 1500 1510 1520 1530 1540 1550 1560 1570 1580 1590 1600 1610 1620 1630 1640 1650 1660 1670 1680 1690 1700 1710 1720 1730 1740 1750 1760 1770 1780 1790 1800 1810 1820 1830 1840 1850 1860 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990', @ischar);
        addParameter(p, 'parallel', '0', @ischar);
        addParameter(p, 'dataroot', '/mnt/sw01-home01/mbmhscc4/scratch/data/Naming_ECoG/avg');
        addParameter(p, 'metadata', 'cvpartition_10fold_leuven_byanimate.mat');
        parse(p, condition, subjects, holdouts, varargin{:});
        fprintf('condition: %s\n', p.Results.condition);
        fprintf('subjects: %s\n', p.Results.subjects);
        fprintf('holdouts: %s\n', p.Results.holdouts);
        fprintf('parallel: %s\n', p.Results.parallel);
    else
        addRequired(p, 'condition');
        addOptional(p, 'subjects', [1:3,5,7:10]);
        addOptional(p, 'holdouts', 1:10);
        addParameter(p, 'windows', 0:10:1990);
        addParameter(p, 'parallel', 0);
        addParameter(p, 'dataroot', 'D:/ECoG/data/avg');
        addParameter(p, 'metadata', 'C:/Users/mbmhscc4/MATLAB/src/MovingWindow_RSA_ECOG/data/cvpartition_10fold_leuven_byanimate.mat');
        parse(p, condition, subjects, holdouts, varargin{:});
        fprintf('condition: %s\n', p.Results.condition);
        fprintf('subjects: %d\n', p.Results.subjects);
        fprintf('holdouts: %d\n', p.Results.holdouts);
        fprintf('parallel: %d\n', p.Results.parallel);
    end

    if isdeployed
        PARALLEL = logical(str2double(p.Results.parallel));
        clean_string = regexprep(p.Results.subjects, ',? *', ' ');
        subjects_set = str2double(strsplit(clean_string));
        clean_string = regexprep(p.Results.holdouts, ',? *', ' ');
        holdout_set = str2double(strsplit(clean_string));
        clean_string = regexprep(p.Results.windows, ',? *', ' ');
        windows_set = str2double(strsplit(clean_string));
        dataroot = p.Results.dataroot;
        metadata_path = p.Results.metadata;
    else
        PARALLEL = p.Results.parallel;
        subjects_set = p.Results.subjects;
        holdout_set = p.Results.holdouts;
        windows_set = p.Results.windows;
        dataroot = p.Results.dataroot;
        metadata_path = p.Results.metadata;
    end

    FIRSTWINDOW = windows_set(1);
    LASTWINDOW = windows_set(end);
    nsubjects = numel(subjects_set);
    nwindows = numel(windows_set);
    nholdout = numel(holdout_set);
    N = nwindows * nsubjects * nholdout;

    [x,y,z] = ndgrid(windows_set,subjects_set,holdout_set);
    [xi,~,~] = ndgrid(1:nwindows,1:nsubjects,1:nholdout);
    windows  = x(:);
    subjects = y(:);
    holdout  = z(:);
    windows_ix  = xi(:);
%     subjects_ix = yi(:);
%     holdout_ix  = zi(:);
    
    if PARALLEL
        ppp = parpool('local');
    end
    % This cv structure was designed to ensure that, when animates and
    % inanimates are modeled at the same time, the holdout sets are composed of
    % a balanced number of animate and inanimate items.
    load(metadata_path, 'CV');

    nchar = 0;
    ps = 0;
    ph = 0;
    
    M = zeros(nwindows,nwindows);
    for i = 1:N
        s = subjects(i);
        w = windows(i);
        h = holdout(i);
        wi = windows_ix(i);
        if w == FIRSTWINDOW;
            resultfile = sprintf('results_%s_%02d_%02d.mat', condition, s, h);
            load(resultfile, 'results', 'fitObj');
        end
        fprintf(repmat('\b',1,nchar));
        nchar = fprintf('% 6.2f%%  s: %02d ps: %02d h: %02d ph: %02d index: %d', 100 * (i / numel(subjects)), s, ps, h, ph, i);
        %if h == 1 && s < 10
        %  continue;
        %end
        for wj = 1:nwindows
            W = windows_set(wj);
            datadir = fullfile(dataroot,'BoxCar','001','WindowStart',sprintf('%04d', W),'WindowSize','0050');
            metafile = fullfile(datadir, 'metadata_raw.mat');
            datafile = fullfile(datadir, sprintf('s%02d_raw.mat', s));
            load(metafile);
            load(datafile);
            
            Mz = selectbyfield(metadata, 'subject', s);
            Fz = selectbyfield(Mz.filters, 'label', 'Leuven');
            leuvenfilter = Fz.filter(:);

            Tz = selectbyfield(Mz.targets, 'label', 'animate');
            switch condition
                case 'animate'
                    leuvenfilter = leuvenfilter & (Tz.target(:) == 1);
                case 'inanimate'
                    leuvenfilter = leuvenfilter & (Tz.target(:) == 0);
                otherwise
                    % pass
            end
            Fz = selectbyfield(Mz.filters, 'label', 'rowfilter');
            rowfilter = Fz.filter(:) & leuvenfilter(:);

            Fz = selectbyfield(Mz.filters, 'label', 'Leuven');
            leuvenfilter = Fz.filter(:);
            cvind = zeros(100, 1);
            for k = 1:size(CV, 2) %#ok<NODEF>
                z = leuvenfilter;
                z(z) = CV(:, k);
                cvind(z) = k;
            end

            Tz = selectbyfield(Mz.targets, ...
                'label'     , 'semantic', ...
                'type'      , 'similarity', ...
                'sim_source', 'Leuven');
            S = Tz.target;
            C = sqrt_truncate_r(S, 0.3);

            Fz = selectbyfield(Mz.filters, 'label', 'colfilter');
            colfilter = Fz.filter(:);

            C = C(rowfilter, :);
            X = X(rowfilter, colfilter);

            testset  = cvind == h;
            testset = testset(rowfilter);

            C_test = C(testset, :);
            X = bsxfun(@minus, X, mean(X));
            X = bsxfun(@rdivide, X, std(X));
            X_test = X(testset, :);

            Cz_test = glmnetPredict(fitObj(wi), X_test);

            M(wi,wj) = norm(C_test - Cz_test,'fro') ./ norm(C_test,'fro');
%             if (wi == wj) && (M(wi,wj) ~= results(wi).err1)
%                 error('Perfomance computation error!')
%             end
        end
        % Windows are the inner loop, so every time we hit the last window
        % it's time to write out a results structure for a given subject
        % and holdout.
        if w == LASTWINDOW;
            resultfile = sprintf('accuracyprofile_%s_%02d_%02d.mat', condition,s,h);
            save(resultfile, 'M', 'windows');
            ps = s;
            ph = h;
        end
    end
    fprintf('\n');

    if PARALLEL
        delete(ppp);
    end
end
function s = selectbyfield(S, varargin)
% SELECTBYFIELD Subset structured array matching on fields within array
% USAGE
% S = struct('label', {'A','A','B','B'}, 'dim', {1,2,1,2});
% s = selectbyfield(S, 'label', 'A');
% s = selectbyfield(S, 'label', 'A', 'dim', 1);
    args = reshape(varargin, 2, []);
    z = true(1, numel(S));
    for i = 1:size(args,2)
        field = args{1, i};
        value = args{2, i};
        if isnumeric(value)
            z = z & cellfun(@(x) isequal(value, x), {S.(field)});
        else
            z = z & strcmp(value, {S.(field)});
        end
    end
    s = S(z);
end
