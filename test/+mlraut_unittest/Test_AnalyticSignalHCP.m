classdef Test_AnalyticSignalHCP < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 11-Nov-2024 23:22:32 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 24.1.0.2689473 (R2024a) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        testObj
    end    

    methods (Static)
        function testObj = called(testObj)
            call_subject(testObj);
        end
    end

    methods (Test)
        function test_afun(this)
            import mlraut.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end

        function test_ctor(this)
            as = this.testObj;
            this.verifyEqual(as.num_nets, 9);
            this.verifyEqual(as.num_sub, 1);
            this.verifyEqual(as.num_tasks, 1);
        end

        function test_task_mask_niigz(this)
            as = this.testObj;
            ic = as.task_mask_niigz;
            % as.task_ref_niigz.view(ic);

            this.verifyEqual(ic.filename, "wmparc.2_binarized.nii.gz")
            this.verifyEqual(ic.qfac, -1)
            this.verifyEqual(size(ic), [91, 109, 91])
            this.verifyInstanceOf(ic.imagingFormat.img, "single")
            this.verifyEqual(dipmax(ic), 1)
            this.verifyEqual(dipsum(ic), 195505)
        end
        
        function test_task_ref_niigz(this)
            as = this.testObj;
            ic = as.task_ref_niigz;
            % ic.view_qc(as.task_mask_niigz);

            this.verifyEqual(ic.filename, "rfMRI_REST1_RL_SBRef.nii.gz")
            this.verifyEqual(ic.qfac, -1)
            this.verifyEqual(size(ic), [91, 109, 91])
            this.verifyInstanceOf(ic.imagingFormat.img, "single")
            this.verifyEqual(dipmax(ic), 49664.171875, AbsTol=1e-3)
            this.verifyEqual(dipsum(ic), 2.286530000098282e+09, AbsTol=1)
        end

        function test_templates(this)
            as = this.testObj;
            as.current_subject = '995174';
            as.current_task = 'rfMRI_REST1_RL';
            
            % template_niigz ~ wmparc
            this.verifyTrue(isfile(as.wmparc_fqfn))
            this.verifyInstanceOf(as.template_niigz, "mlfourd.ImagingContext2")
            this.verifyEqual(as.template_niigz.filename, "wmparc.2.nii.gz")

            % template_cifti ~ task_ref_dscalar_fqfn
            this.verifyTrue(isfile(as.task_ref_dscalar_fqfn));
            this.verifyTrue(isstruct(as.template_cifti.metadata));
            this.verifyTrue(iscell(as.template_cifti.diminfo));
            this.verifyEqual(as.template_cifti.diminfo{1}.models{1}.count, 29696);
            this.verifyEqual(as.template_cifti.diminfo{1}.models{1}.struct, 'CORTEX_LEFT');
            this.verifyEqual(as.template_cifti.diminfo{1}.models{1}.type, 'surf');
            this.verifyEqual(size(as.template_cifti.diminfo{1}.models{1}.vertlist), [1, 29696]);
            this.verifyEqual(as.template_cifti.diminfo{1}.models{2}.count, 29716);
            this.verifyEqual(as.template_cifti.diminfo{1}.models{2}.struct, 'CORTEX_RIGHT');
            this.verifyEqual(as.template_cifti.diminfo{1}.models{2}.type, 'surf');
            this.verifyEqual(size(as.template_cifti.diminfo{1}.models{2}.vertlist), [1, 29716]);
            this.verifyEqual(size(as.template_cifti.cdata), [91282, 1]);
        end

        function test_average_network_signals(this)
            %% candidate supplemental figure, with concat runs

            as = this.testObj;
            call_subject(as);

            yeos = as.HCP_signals;
            yeo_names = mlraut.NetworkData.NETWORKS_YEO_NAMES;
            for anat = ["cbm", "ctx", "str", "thal"]
                figure
                hold on
                for y = 1:length(yeo_names)
                    plot(real(yeos.(anat).psi(:, y))); 
                end
                title(sprintf("yeos.%s", anat));  
                legend(yeo_names)
                hold off

                figure;     
                tiledlayout(3,3);
                for y = 1:length(yeo_names)
                    nexttile
                    histogram(real(yeos.(anat).psi(:, y))); 
                    title(sprintf("yeos(%s)", yeo_names{y}), 100);
                    title(sprintf("yeos.%s(%s)", anat, yeo_names{y}));    
                end
            end
        end
        






        function test_fultz_iFV(this)
            as = this.testObj;
            as.source_physio = "iFV";
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "Y", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_fultz_sFV(this)
            as = this.testObj;
            as.source_physio = "sFV";
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "Y", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_fultz_latV(this)
            as = this.testObj;
            as.source_physio = "latV";
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "Y", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_fultz_csf(this)
            as = this.testObj;
            as.source_physio = "csf";
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "Y", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end



        

        function test_call_wmparc(this)

            return

            % wmparc = 'precuneus';
            % wmparc = 'posteriorcingulate';
            % wmparc = 'hippocampus';
            % wmparc = 'entorhinal';
            % wmparc = 'medialorbitofrontal';
            % wmparc = 'insula';
            % wmparc = 'cuneus';
            % wmparc = 'thalamus';
            % wmparc = 'caudate';
            % wmparc = 'putamen';
            % wmparc = 'pallidum';
            % wmparc = 'cerebellum';
            % wmparc = 'ponsvermis';
            % wmparc = 'brainstem';
            % wmparc = 'brainstem+';
            % wmparc = 'csf';
            % wmparc = 'centrumsemiovale';
            % wmparc = 'corpuscallosum';

            wmparcs = { ...
                'cuneus' 'corpuscallosum' ...
                'precuneus' 'posteriorcingulate' 'hippocampus' 'entorhinal' 'medialorbitofrontal' ...
                'insula' ...
                'thalamus' 'caudate' 'putamen' 'pallidum' 'cerebellum' ...
                'ponsvermis' 'brainstem' 'brainstem+' 'csf' 'centrumsemiovale'};

            for w = wmparcs
                wmparc = w{1};

                as = this.testObj;
                as.do_save=false;
                as.do_save_dynamic=false;
                as.do_save_ciftis=true;
                as.source_physio=wmparc;
                as.out_dir = sprintf('/Volumes/PrecunealSSD2/AnalyticSignalHCP/physio_%s', wmparc);

                disp(as)
                call(as);
            end
        end
    end
    
    methods (TestClassSetup)
        function setupAnalyticSignalHCP(this)
            import mlraut.*
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            this.testObj_ = AnalyticSignalHCP( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_global_signal_regression=true, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                v_physio=50, ...
                plot_range=1:250, ...
                source_physio="iFV", ...
                tags=stackstr(use_dashes=true));
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalHCPTest(this)
            this.testObj = copy(this.testObj_);
            malloc(this.testObj);
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
