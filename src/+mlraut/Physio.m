classdef Physio < handle 
    %% Implements Raut, et al.  Global waves synchronize the brain's functional systems with fluctuating arousal.
    %  https://www.science.org/doi/10.1126/sciadv.abf2709
    %  Extends physio_phase_mapping.m
    %  
    %  Created 25-Jul-2022 11:54:11 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.12.0.1975300 (R2022a) Update 3 for MACI64.  Copyright 2022 John J. Lee.
  
    properties
        do_save
        source_physio
        gs_subtract % false -> normalize by global signal

        hp_thresh = 0.01; % lower bound, Ryan ~ 0.01, converted to units of 1/tr
        lp_thresh = 0.1; % higher bound, Ryan ~ 0.05-0.1, CSF ~ 0.1, converted to units of 1/tr

        abs_h1 % amplitude of arousal
        abs_h2 % amplitude of phase-locked analytic signal
        angles_m2pi % phase mod 2pi
        bold
        h1
        h2
        physio_ds % physiological signal resampled to match BOLD
        plvs_xt % phase-locked BOLD in spacetime

        current_subject = ''
        current_task = ''
    end

    properties (Dependent)
        dmn_parcels
        num_sub
        num_tasks
        subjects
        tasks    
    end

    methods %% GET
        function g = get.dmn_parcels(~)
            g = {'L_S_parieto_occipital', 'L_G_precuneus', 'L_S_subparietal', ...
                 'R_S_parieto_occipital', 'R_G_precuneus', 'R_S_subparietal'};
        end
        function g = get.num_sub(this)
            g = numel(this.subjects);
        end
        function g = get.num_tasks(this)
            g = numel(this.tasks);
        end
        function g = get.subjects(~)
            if contains(hostname, 'precuneal')
                g = {'995174'};
                return
            end
            if contains(hostname, 'pascal')
                g = {'995174'};
                return
            end
            if contains(hostname, 'cluster')
                g = {'995174'};
                return
            end
            error('mlraut:NotImplementedError', 'Physio.get.subjects');
        end
        function g = get.tasks(this)
            g = this.tasks_;
        end
    end

    methods
        function this = average_network_signals(this, BOLD, s, t)
            assert(isscalar(s));
            assert(isscalar(t));

            for n = 1:7
                this.ctx_signals(:,n,s,t) = mean(BOLD(:,this.mask_ctx_HCP & this.networks_HCP==n),2,'omitnan'); %#ok<*mean> 
                this.str_signals(:,n,s,t) = mean(BOLD(:,this.mask_str_HCP & this.networks_HCP==n),2,'omitnan');
                this.thal_signals(:,n,s,t) = mean(BOLD(:,this.mask_thal_HCP & this.networks_HCP==n),2,'omitnan');
                this.cbm_signals(:,n,s,t) = mean(BOLD(:,this.mask_cbm_HCP & this.networks_HCP==n),2,'omitnan');
            end
            this.ctx_signals(:,8,s,t) = mean(BOLD(:,this.mask_ctx_HCP),2,'omitnan');
            this.str_signals(:,8,s,t) = mean(BOLD(:,this.mask_str_HCP),2,'omitnan');
            this.thal_signals(:,8,s,t) = mean(BOLD(:,this.mask_thal_HCP),2,'omitnan');
            this.cbm_signals(:,8,s,t) = mean(BOLD(:,this.mask_cbm_HCP),2,'omitnan');
        end
        function this = call(this)

            for s = 1:this.num_sub
                tic
                        
                disp(['Processing subject ' this.subjects{s} ' : ' num2str(s) ' out of ' num2str(this.num_sub)]);
                subj = num2str(this.subjects{s});            
                
                for t = 1:this.num_tasks
                    task = this.tasks{t};            
                    disp(['Processing task ' task ' : ' num2str(t) ' out of ' num2str(this.num_tasks)]);

                    % BOLD
                    try
                        BOLD = this.task_dtseries(subj, task); 
                        assert(~isempty(BOLD))
                    catch ME
                        disp([subj ' ' task ' BOLD missing or defective:']);
                        disp(ME)
                        continue
                    end
                     
                    % Physio
                    try
                        physio = this.task_physio(subj, task, BOLD);
                    catch ME
                        disp([subj ' ' task ' physio missing:']);
                        disp(ME)
                        continue
                    end
            
                    %if size(BOLD,1)~=1191, continue, end
            
                    % Center and rescale
                    physio = this.center_and_rescale(physio);
                    BOLD = this.global_signal_regression(BOLD);
                    BOLD = this.center_and_rescale(BOLD);

                    % Filtering
                    if ~isempty(this.lp_thresh) && ~isempty(this.hp_thresh)
                        fprintf('Filter physio_ds & BOLD by: [%g, %g]\n', this.hp_thresh, this.lp_thresh)
                        [b,a] = butter(2,[this.hp_thresh,this.lp_thresh]/(this.Fs/2));
                        if ~isempty(physio)
                            physio = single(filtfilt(b,a,double(physio)));
                        end
                        BOLD = single(filtfilt(b,a,double(BOLD)));
                    end
            
                    % Phase-locking values
                    if ~isempty(physio)
                        h1_ = hilbert(physio);
                        this.h1(:,s,t) = h1_;
                        this.abs_h1(:,s,t) = abs(h1_);
                        this.physio_ds(:,s,t) = physio;
                    end
                    h2_ = hilbert(BOLD);
                    this.bold(:,:,s,t) = BOLD;
                    this.h2(:,:,s,t) = h2_;
                    this.abs_h2(:,:,s,t) = abs(h2_);
                    if isempty(physio)
                        angles_ = -unwrap(angle(h2_));
                        this.angles_m2pi(:,:,s,t) = mod(angles_, 2*pi);
                        this.plvs_xt(:,:,s,t) = real(h2_);
                    else
                        angles_ = bsxfun(@minus,unwrap(angle(h2_)),unwrap(angle(h1_))); % unwraps to -inf:inf
                        this.angles_m2pi(:,:,s,t) = mod(angles_, 2*pi);
                        this.plvs_xt(:,:,s,t) = real(h2_./h1_);
                    end
                    plvs_ = mean(this.plvs_xt(:,:,s,t),'omitnan');

                    this.write_cifti(this.bold(:,:,s,t), sprintf('bold_%i_%i', s, t));
                    this.write_cifti(this.abs_h2(:,:,s,t), sprintf('abs_h2_%i_%i', s, t));
                    this.write_cifti(this.angles_m2pi(:,:,s,t), sprintf('angles_m2pi_%i_%i', s, t));
                    this.write_cifti(this.plvs_xt(:,:,s,t), sprintf('plvs_xt_%i_%i', s, t));
                    this.write_cifti(plvs_, sprintf('plvs_%i_%i', s, t));                            

                    % Get mean network signals
                    this.average_network_signals(this.plvs_xt(:,:,s,t), s, t);

                    toc                    
                end
            end

            if this.do_save
                save(fullfile(this.out_dir, 'mlraut_Physio_call.mat'), 'this');
            end
        end
        function dat = center_and_rescale(~, dat)
            if isempty(dat)
                return
            end
            dat = dat - mean(dat, 'all', 'omitnan');
            dat = dat / max(abs(dat), [], 'all');
        end
        function dat = global_signal_regression(this, dat)
            assert(~isempty(dat))
            assert(size(dat, 1) == this.num_frames)
            assert(size(dat, 2) == this.num_nodes)

            if this.gs_subtract
                dat = dat - mean(dat, 2);
            else
                dat = dat ./ mean(dat, 2);
            end
        end
        function physio = task_physio(this, subj, task, BOLD)
            switch char(this.source_physio)
                case 'RV'
                    physio = this.physio_rv(subj, task, BOLD);
                    physio = this.trim_frames(physio);
                case 'HRV'
                    physio = this.physio_hrv(subj, task, BOLD);
                    physio = this.trim_frames(physio);
                case {'iFV', 'iFV-brightest'}
                    physio = this.physio_iFV(subj, task);
                    physio = this.trim_frames(physio);
                otherwise
                    physio = [];
            end
        end
        function physio = physio_hrv(this, subj, task, BOLD)
        end
        function physio = physio_iFV(this, subj, task)
            this.current_subject = subj;
            this.current_task = task;
            iFV = mlraut.IFourthVentricle(this);
            physio = call(iFV);
        end
        function physio = physio_rv(this, subj, task, BOLD)
            %  Args:
            %      subj (text)
            %      task (text)
            %      BOLD (numeric)
            %  Returns:
            %      physio (numeric): from <task>_Physio_log.txt

            assert(istext(subj));
            assert(istext(task));
            assert(isnumeric(BOLD));
            fqfn = fullfile( ...
                this.data_dir(subj, task), strcat(task, '_Physio_log.txt'));

            data_ = importdata(fqfn);

            % Get 6 sec windows
            assert(length(data_)/400 >= 860, 'task_physio().data_')

            time_vec_bold = this.tr*(1:size(BOLD,1))';
            time_vec_phys = (0:length(data_)-1)'/400;
            physio = zeros(size(time_vec_bold));
            
            data_(:,3) = zscore(data_(:,3));
            for i = 5:length(physio)-4
                
                % For RV
                [~,phys_start] = min(abs(time_vec_phys-(time_vec_bold(i)-3)));
                [~,phys_end] = min(abs(time_vec_phys-(time_vec_bold(i)+3)));
                physio(i) = std(data_(phys_start:phys_end,2));
               
                % For HRV
                %[pks,locs] = findpeaks(physio(phys_start:phys_end,3),'minpeakdistance',round(400/(180/60)));
                %             %,'minpeakwidth',400/(1/(200/60))); % max heart rate = 180 bpm; at 400 Hz, minimum of 100 samples apart
                %locs = locs(pks>prctile(physio(phys_start:phys_end,3),60));
                %Avg1(i) = mean(diff(locs),'omitnan')/400;
            end
        end
        
        function this = Physio(varargin)
            %% PHYSIO 
            %  Args:
            %      out_dir (folder): for physiological data intermediates, default is pwd.
            %      hp_thresh (isnumeric): default := 0.01 from Ryan.
            %      lp_thresh (isnumeric): default := 0.05-0.1 from Ryan.
            %      do_save (logical): save mlraut.Pysio to mat ~ 10 GB.
            %      tasks (cell of text): '', 'iFV-brightest', 'iFV', 'RV', 'HRV'
            %      source_physio (logical)
            %      gs_subtract (logical)

            error("mlraut:DeprecationError", stackstr())
            
            %this.tasks_ = {'rfMRI_REST1_7T_AP','rfMRI_REST1_7T_PA' 'rfMRI_REST2_7T_AP','rfMRI_REST2_7T_PA'};
            this.tasks_ = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
            %this.tasks_ = {'tfMRI_MOTOR_LR','tfMRI_MOTOR_RL'};

            ip = inputParser;
            addParameter(ip, "out_dir", pwd, @isfolder);
            addParameter(ip, "hp_thresh", this.hp_thresh, @isnumeric);
            addParameter(ip, "lp_thresh", this.lp_thresh, @isnumeric);
            addParameter(ip, "do_save", true, @islogical);
            addParameter(ip, "tasks", this.tasks_, @iscell);
            addParameter(ip, "source_physio", 'iFV-brightest', @istext)
            addParameter(ip, "gs_subtract", true, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.out_dir_ = ipr.out_dir;
            this.hp_thresh = ipr.hp_thresh;
            this.lp_thresh = ipr.lp_thresh;
            this.do_save = ipr.do_save;
            this.tasks_ = ipr.tasks;
            this.source_physio = ipr.source_physio;
            this.gs_subtract = ipr.gs_subtract;

            this.ctx_signals = single(nan(this.num_frames,8,this.num_sub,this.num_tasks));
%            this.prec_signals = single(nan(this.num_frames,6,this.num_sub,this.num_tasks));
            this.str_signals = single(nan(this.num_frames,8,this.num_sub,this.num_tasks));
            this.thal_signals = single(nan(this.num_frames,8,this.num_sub,this.num_tasks));
            this.cbm_signals = single(nan(this.num_frames,8,this.num_sub,this.num_tasks));

            this.physio_ds = single(nan(this.num_frames, this.num_sub, this.num_tasks));
            this.bold = single(nan(this.num_frames, this.num_nodes, this.num_sub, this.num_tasks));
            this.h1 = single(nan(this.num_frames, this.num_sub, this.num_tasks));
            this.h2 = single(nan(this.num_frames, this.num_nodes, this.num_sub, this.num_tasks));
            this.abs_h1 = single(nan(this.num_frames, this.num_sub, this.num_tasks));
            this.abs_h2 = single(nan(this.num_frames, this.num_nodes, this.num_sub, this.num_tasks));
            this.plvs_xt = single(nan(this.num_frames, this.num_nodes, this.num_sub, this.num_tasks));
            this.angles_m2pi = single(nan(this.num_frames, this.num_nodes, this.num_sub, this.num_tasks));

            addpath(genpath(fullfile(this.waves_dir, 'Dependencies', '')));
            addpath(genpath(fullfile(this.waves_dir, 'supporting_files', '')));
        end
    end

    methods (Static)
        function this = sweep_spectral_range(varargin)
            %  Params:
            %      tag (text): folders are named "arousal-waves-<tag>-0p01-0p05"
            %      physio (text): '', 'iFV-brightest', 'iFV', 'HRV', 'RV'
            %      gs_subtract (logical)
            %      fmin (scalar): > 0, units of 1/tr.
            %      fmax (scalar): < inf, units of 1/tr.
            %      N (optional scalar): num. of requested samples of f, default is 9.

            ip = inputParser;
            addParameter(ip, "tag", "iFV-brightest", @istext);
            addParameter(ip, "physio", "iFV-brightest", @istext);
            addParameter(ip, "gs_subtract", true, @islogical)
            addParameter(ip, "fmin", 0, @isscalar);
            addParameter(ip, "fmax", inf, @isscalar);
            addParameter(ip, "N", 9, @isscalar);
            parse(ip, varargin{:});
            ipr = ip.Results;
            if ipr.fmin < 2/1191
                ipr.fmin = 1/512;
            end
            if ipr.fmax > 0.5
                ipr.fmax = 0.5;
            end

            % estimate df ~ log sweep of spectral range
            Dlogf = log2(ipr.fmin) - log2(ipr.fmax);
            N = ipr.N - 1;
            set_f = [flip(2.^(log2(ipr.fmax):Dlogf/N:log2(ipr.fmin)))];

            disp(set_f')
            prefix = 'arousal-waves';
            if ~ipr.gs_subtract
                prefix = strcat(prefix, '-gsq'); % global signal quotient
            end
            if "" ~= ipr.tag
                prefix = strcat(prefix, '-', ipr.tag);
            end

            % wide-band or no filtering
            fold = strcat(prefix, '-0.002-0.5');
            if ~isfolder(fold)
                mkdir(fold);
                pwd0 = pushd(fold);
                fprintf('mlraut.Physio.sweep_spectral_range:  working in %s\n', pwd)
                this = mlraut.Physio('hp_thresh', 1/512, 'lp_thresh', 0.5, ...
                    'source_physio', ipr.physio, ...
                    'gs_subtract', ipr.gs_subtract);
                call(this)
                popd(pwd0);
            end

            % sweep filter bands in calls to Physio
            for idx = 1:length(set_f)-1
                fold = sprintf('%s-%g-%g', ...
                    prefix, round(set_f(idx),2,'significant'), round(set_f(idx+1),2,'significant'));
                fold = strrep(fold, '.', 'p');
                if ~isfolder(fold)
                    mkdir(fold);
                    pwd0 = pushd(fold);
                    fprintf('mlraut.Physio.sweep_spectral_range:  working in %s\n', pwd)
                    this = mlraut.Physio('hp_thresh', set_f(idx), 'lp_thresh', set_f(idx+1), ...
                        'source_physio', ipr.physio, ...
                        'gs_subtract', ipr.gs_subtract);
                    call(this)
                    popd(pwd0);
                end
            end
        end
    end

    %% PRIVATE

    properties (Access = private)
        subjects_
        tasks_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
