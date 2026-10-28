close all;
study_dir=fileparts(mfilename('fullpath'));
addpath(study_dir,'-begin');
if ~exist('ablation_options','var'), ablation_options=struct(); end

T=1/50; duration=30;
if isfield(ablation_options,'duration'), duration=ablation_options.duration; end
assert(isscalar(duration) && duration>=T && isfinite(duration));
cfg=struct();
cfg.version='figure-eight-common-baselines-20261019';
cfg.dt=T; cfg.t=0:T:duration; cfg.scenarioId=3;
cfg.seed=21261026;
if isfield(ablation_options,'seed'), cfg.seed=ablation_options.seed; end
cfg.g=[0;0;-9.81]; cfg.mUav=1.5059; cfg.mPayload=0.1382;
cfg.cableLength=0.500;
cfg.sensor=struct('bearingSigma',deg2rad(0.02568956),'bearingBias', ...
    [3.41866553510e-3;6.53223936837e-4;-2.55318581495e-4], ...
    'bearingNoiseShape',eye(3), ...
    'bearingCorrelation',0.80,'visualEnvelopeMinimum',0.85, ...
    'visualEnvelopeAmplitude',0.30,'uavPositionSigma',0.004, ...
    'cableSigma',0.004,'imuSigma',0.005,'thrustRelativeSigma',0.001, ...
    'attitudeSigma',deg2rad(0.08),'outlierProbability',0.0027,'outlierScale',13);
cfg.nominal=struct('initialPositionSigma',0.004,'initialBearingSigma',deg2rad(1.5), ...
    'uavPositionVariance',6.10e-5,'cableVariance',6.10e-5, ...
    'angularVariance',3.59e-3);
cfg.shf=struct('rho',0.92);
cfg.initialRange=0.35; cfg.initialRangeSigma=0.25;
cfg.initialVelocitySigma=0.20; cfg.modelAccelerationSigma=0.80;
cfg.vb=struct('iterations',2,'rho',0.90,'alpha',6, ...
    'equivalentSamples',30,'qScale',0.3,'rScale',0.1);
if isfield(ablation_options,'selectedParameters')
    hp=ablation_options.selectedParameters;
    cfg.vb=struct('iterations',hp.Iterations,'rho',hp.Rho,'alpha',hp.Alpha, ...
        'equivalentSamples',hp.EquivalentSamples,'qScale',hp.QScale,'rScale',hp.RScale);
end
assert(cfg.vb.iterations>=1 && cfg.vb.iterations==floor(cfg.vb.iterations));
assert(cfg.vb.rho>0 && cfg.vb.rho<=1 && cfg.vb.alpha>0 && ...
    cfg.vb.equivalentSamples>0 && cfg.vb.qScale>0 && cfg.vb.rScale>0);
cfg.epsilon=1e-12;
cfg.timeAverageEnvelopeSquared=mean((0.5+2.5*sin(2*pi*0.12*cfg.t+0.4).^2).^2);
cfg.timeAverageVisualEnvelopeSquared=mean((cfg.sensor.visualEnvelopeMinimum+ ...
    cfg.sensor.visualEnvelopeAmplitude*sin(2*pi*0.12*cfg.t+0.4).^2).^2);

cfg.lengthBias_m=[0.05,0.10];
cfg.expectedEpsilonC_m=-cfg.lengthBias_m;
cfg.bisgaardLengthErrorSD_m=0.05/sqrt(3);
cfg.bisgaardLengthSeedOffset=420000;
make_figures=true; save_figures=true; save_results=true;
if isfield(ablation_options,'makeFigures'), make_figures=ablation_options.makeFigures; end
if isfield(ablation_options,'saveFigures'), save_figures=ablation_options.saveFigures; end
if isfield(ablation_options,'saveResults'), save_results=ablation_options.saveResults; end
out_dir=fullfile(study_dir,'Results');
if isfield(ablation_options,'outputDir'), out_dir=ablation_options.outputDir; end
if save_results || (make_figures && save_figures)
    if ~isfolder(out_dir), mkdir(out_dir); end
end

Variant=["VBF-CL";"without PAE";"without CL";"without VB"; ...
    "VBF-CL, E[epsilon_c(k)] = -0.05 m";"VBF-CL, E[epsilon_c(k)] = -0.10 m"; ...
    "Bisgaard-UKF";"PCMPC-EKF"];
plot_labels={'VBF-CL','without PAE','Bearing-only','without VB', ...
    'VBF-CL, $\mathrm{I\!E}[\epsilon_c(k)]=-0.05\,\mathrm{m}$', ...
    'VBF-CL, $\mathrm{I\!E}[\epsilon_c(k)]=-0.10\,\mathrm{m}$', ...
    'Bisgaard-UKF','PCMPC-EKF'};
UsePAE=[true;false;true;true;true;true;false;false];
UseCL=[true;true;false;true;true;true;true;true];
UseVB=[true;true;true;false;true;true;false;false];
LengthBias_m=[0;0;NaN;0;cfg.lengthBias_m(:);0;0];
ExpectedEpsilonC_m=-LengthBias_m;
MeasurementModel=repmat("CL-pseudolinear",8,1);
MeasurementModel(3)="bearing-only";
MeasurementModel(7)="pendulum bearing UKF";
MeasurementModel(8)="body-frame cable position/rate EKF";
InitialPrior=[repmat("common six-state prior",6,1); ...
    "Experiment_4 UKF prior";"Experiment_4 EKF prior"];
variant_table=table(Variant,UsePAE,UseCL,UseVB,LengthBias_m,ExpectedEpsilonC_m,MeasurementModel,InitialPrior);
variant_table.ExtraLengthErrorSD_m=[0;0;NaN;0;0;0;cfg.bisgaardLengthErrorSD_m;0];
n_methods=height(variant_table);
Additional_Variant=["EKF";"EKF without PAE";"PL-KF";"PL-KF without PAE"; ...
    "SHF";"SHF without PAE";"Bearing-only"];
additional_functions={@EKF,@EKF_without_PAE,@PLKF,@PLKF_without_PAE,@SHF,@SHF_without_PAE,@Bearing_only};
Additional_status=cell(1,numel(Additional_Variant));
Additional_parameters=cell(1,numel(Additional_Variant));
algorithm_parameter_table=table(Variant,UseVB,repmat(cfg.vb.iterations,n_methods,1), ...
    repmat(cfg.vb.rho,n_methods,1),repmat(cfg.vb.alpha,n_methods,1), ...
    repmat(cfg.vb.equivalentSamples,n_methods,1),repmat(cfg.vb.qScale,n_methods,1), ...
    repmat(cfg.vb.rScale,n_methods,1),'VariableNames', ...
    {'Variant','UseVB','Iterations','Rho','Alpha','n0','c_Q','c_R'});
algorithm_parameter_table{~UseVB,3:6}=NaN;
algorithm_parameter_table{7:8,7:8}=NaN;
example=struct();
run_clock=tic;
run_seed=cfg.seed+20000*cfg.scenarioId+1;
rng(run_seed,'twister');
data=generate_simulation(cfg);
observations=data.observations;
result=cell(n_methods,1);
result{1}=VBF_CL(observations,cfg);
result{2}=VBF_without_PAE(observations,cfg);
result{3}=VBF_without_CL(observations,cfg);
result{4}=VBF_without_VB(observations,cfg);
for bias_id=1:numel(cfg.lengthBias_m)
    result{4+bias_id}=VBF_CL_uncertain_length(observations,cfg,cfg.lengthBias_m(bias_id));
end
estimation_idx=1:numel(cfg.t);
[l_Bisgaard,Bisgaard_length]=Cable_length_error( ...
    cfg.cableLength,numel(cfg.t),run_seed+cfg.bisgaardLengthSeedOffset, ...
    cfg.bisgaardLengthErrorSD_m,estimation_idx);
[X_Bisgaard_UKF,Bisgaard_parameters]=Bisgaard_UKF( ...
    observations.bearing,observations.pUav,observations.pUav, ...
    observations.accelerationUav,T,l_Bisgaard,estimation_idx);
[X_PCMPC_EKF,PCMPC_parameters]=PCMPC_EKF( ...
    observations.bearing,observations.pUav,observations.pUav, ...
    observations.attitude,observations.thrust,T,cfg.cableLength, ...
    cfg.mUav,estimation_idx);
result{7}=struct('X',X_Bisgaard_UKF);
result{8}=struct('X',X_PCMPC_EKF);
additional_result=cell(numel(Additional_Variant),1);
for method=1:numel(Additional_Variant)
    additional_result{method}=additional_functions{method}(observations,cfg);
    br=additional_result{method};
    assert(isequal(br.initialState,result{1}.initialState) && ...
        isequal(br.initialCovariance,result{1}.initialCovariance) && ...
        isequal(br.Q,result{1}.Q),'Common initial state/P0/Q differ.');
    Additional_parameters{1,method}=rmfield(br,'X');
    if isfield(br,'Failure'),Additional_status{1,method}=br.Failure;
    else,Additional_status{1,method}='';end
end
for pair=[1 3 5]
    one=additional_result{pair};two=additional_result{pair+1};
    assert(isequal(one.Rinitial,two.Rinitial),'PAE pair Rinitial differs.');
end
assert(isequaln(additional_result{3}.X,result{4}.X),'PL-KF must match without VB.');
for method=1:n_methods
    if method<=6
        assert(isequal(result{method}.initialState,result{1}.initialState) && ...
            isequal(result{method}.initialCovariance,result{1}.initialCovariance), ...
            'Initial priors differ across the six ablations.');
    end
end

example.data=data;
example.X=zeros(6,numel(cfg.t),n_methods);
for method=1:n_methods, example.X(:,:,method)=result{method}.X; end
example.seed=run_seed;
example.AdditionalX=nan(6,numel(cfg.t),numel(Additional_Variant));
for method=1:numel(Additional_Variant)
    example.AdditionalX(:,:,method)=additional_result{method}.X;
end
Common_parameters=cell(11,4);
common_result=[result(1:4);additional_result];
for method=1:11
    br=common_result{method};
    Common_parameters(method,:)={br.initialState,br.initialCovariance,br.Q,br.Rinitial};
end
Common_parameter_table=table([Variant(1:4);Additional_Variant], ...
    Common_parameters(:,1),Common_parameters(:,2),Common_parameters(:,3),Common_parameters(:,4), ...
    'VariableNames',{'Variant','InitialState','P0','Q','Rinitial'});

literature_parameter_table=table(Variant(7:8), ...
    {Bisgaard_parameters.P0;PCMPC_parameters.P0}, ...
    {Bisgaard_parameters.Q;PCMPC_parameters.Q}, ...
    {Bisgaard_parameters.R;PCMPC_parameters.R}, ...
    [Bisgaard_parameters.derivative_tau;PCMPC_parameters.derivative_tau], ...
    'VariableNames',{'Variant','P0','Q','R','DerivativeTau_s'});
colors=[0 .18 .65; .90 .40 0; 0 .60 .20; .55 .20 .70; .85 .55 .05; .8 .15 .15; ...
    0 .55 .75; .75 .20 .55];
styles=repmat({'-'},1,n_methods);
figure_method_order=[2:n_methods,1];
shown=example; t=shown.data.t;
letters='abcdefghijklmnop';
Figure_Index=table(5,"Eight trajectories versus Ground Truth (2x4)", ...
    "Figure_05.fig","ablation_trajectories.pdf", ...
    'VariableNames',{'Figure','Content','EditableFile','ManuscriptFile'});
if make_figures

    figure(5); clf;
    set(gcf,'Color','w','Name','Ablation trajectories','Units','inches', ...
        'Position',[1 1 10.8 4.5]);
    tl=tiledlayout(2,4,'TileSpacing','compact','Padding','loose');
    all_points=[shown.data.pPayloadTrue,reshape(shown.X(1:3,:,:),3,[])];
    limits=[min(all_points,[],2),max(all_points,[],2)];
    margin=.08*(limits(:,2)-limits(:,1));
    limits=limits+[-margin,margin];
    for panel=1:n_methods
        method=figure_method_order(panel);
        ax=nexttile(tl); hold(ax,'on');
        p=shown.data.pPayloadTrue; x=shown.X(1:3,:,method);
        hgt=plot3(p(1,:),p(2,:),p(3,:),'k-','LineWidth',1.0);
        he=plot3(x(1,:),x(2,:),x(3,:),'Color',colors(method,:), ...
            'LineStyle','--','LineWidth',1.05);
        axis equal; xlim(limits(1,:)); ylim(limits(2,:)); zlim(limits(3,:));
        grid on; view(45,25); set(ax,'FontName','Times New Roman','FontSize',8);
        xlabel('X (m)','rotation',-15); ylabel('Y (m)','rotation',25); zlabel('Z (m)');
        if method<=4 || method>=7
            panel_title=sprintf('(%c) %s',letters(panel),plot_labels{method});
        else
            panel_title={sprintf('(%c) VBF-CL',letters(panel)), ...
                sprintf('$\\mathrm{I\\!E}[\\epsilon_c(k)]=%.2f\\,\\mathrm{m}$', ...
                cfg.expectedEpsilonC_m(method-4))};
        end
        title(panel_title,'Interpreter','latex','FontSize',9,'FontWeight','normal');
        if panel==1
            trajectory_legend_ax=ax;
            trajectory_ground_truth_handle=hgt;
        end
    end

    trajectory_legend_handles=gobjects(1,n_methods+1);
    trajectory_legend_handles(1)=trajectory_ground_truth_handle;
    for panel=1:n_methods
        method=figure_method_order(panel);
        trajectory_legend_handles(panel+1)=plot3(trajectory_legend_ax,NaN,NaN,NaN, ...
            'Color',colors(method,:),'LineStyle','--','LineWidth',1.05);
    end
    trajectory_legend=legend(trajectory_legend_ax,trajectory_legend_handles, ...
        [{'Ground Truth'},plot_labels(figure_method_order)],'Interpreter','latex', ...
        'FontName','Times New Roman','FontSize',9,'Orientation','horizontal', ...
        'NumColumns',5,'Box','off','AutoUpdate','off');
    trajectory_legend.ItemTokenSize=[16 10];
    trajectory_legend.Layout.Tile='south';

    if save_figures
        for fig_id=1:height(Figure_Index)
            fig=figure(Figure_Index.Figure(fig_id));
            drawnow;
            savefig(fig,fullfile(out_dir,Figure_Index.EditableFile(fig_id)));
            export_trajectory_figure(fig,fullfile(out_dir,Figure_Index.ManuscriptFile(fig_id)));
        end
    end
end

environment=struct('MATLABVersion',version,'Computer',computer, ...
    'RunTime',char(datetime('now','Format','yyyyMMdd''T''HHmmss')), ...
    'SourceFile',mfilename('fullpath'),'NumRuns',1,'ElapsedSeconds',toc(run_clock));
if save_results
    save(fullfile(out_dir,'Single_Run.mat'),'cfg','example','variant_table', ...
        'algorithm_parameter_table','literature_parameter_table', ...
        'Bisgaard_parameters','PCMPC_parameters','Bisgaard_length', ...
        'Additional_Variant','Additional_parameters','Additional_status', ...
        'Common_parameter_table','Figure_Index','environment');
end

function export_trajectory_figure(fig,path)
warning_state=warning('query','MATLAB:print:ContentTypeImageSuggested');
warning_cleanup=onCleanup(@() warning(warning_state));
warning('off','MATLAB:print:ContentTypeImageSuggested');
exportgraphics(fig,path,'ContentType','vector','BackgroundColor','white');
end
