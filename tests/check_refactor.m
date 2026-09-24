function report = check_refactor(action)
%CHECK_REFACTOR Save/compare numerical results and benchmark representative work.
% The baseline is local test data, not another copy of the simulation file.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
config = orb5_config(struct('path',root,'visible','off'));
baseline = fullfile(root,'tests','refactor_baseline.mat');
oldVisible = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
cleanup = onCleanup(@()finish(oldVisible));

tic;
potential = orb5_plot_library.load_data(config);
report.potential_load_seconds = toc;
tic;
transport = orb5_plot_library.load_flux_data(config);
report.transport_load_seconds = toc;

% Keep derived fields in full, and sample the large raw spectral tensor.
fields = {'t','t_norm','s','nlist','phi_Z','EZ','EZZ','EZ_max', ...
    'radial_envelope','lfs_ES','lfs_ES_max','efield_no_zf','efield_zf'};
snapshot.potential = struct();
for k = 1:numel(fields)
    snapshot.potential.(fields{k}) = potential.(fields{k});
end
snapshot.potential.spectral_sample = potential.p3d(1:100:end,:,:,:);
snapshot.transport = transport;

times = zeros(1,3);
for k = 1:3
    tic;
    entry = orb5_plot_library.find_plot('cross_correlation',config);
    entry.callback(potential,config,round(potential.ns/2));
    times(k) = toc;
    close all force;
end
report.correlation_seconds = median(times);
signal = reshape(potential.lfs_ES(:,end,:),potential.nt,potential.ns);
times = zeros(1,3);
for k = 1:3
    tic;
    [coefficients,frequencies] = cwt_complex(signal,1/potential.dt_norm,config.wname,800);
    times(k) = toc;
end
report.spectrum_seconds = median(times);
snapshot.spectrum = coefficients;
snapshot.frequencies = frequencies;

if strcmp(action,'save')
    save(baseline,'snapshot','report','-v7.3');
    fprintf('Saved numerical baseline.\n');
elseif strcmp(action,'compare')
    previous = load(baseline);
    compare(previous.snapshot,snapshot,'snapshot');
    disp('PASS: refactored data and spectrum match the pre-refactor baseline.');
    disp('Before:'); disp(previous.report);
    disp('After:'); disp(report);
    fid = fopen(fullfile(root,'tests','refactor_benchmark.json'),'w');
    fprintf(fid,'%s',jsonencode(struct('before',previous.report,'after',report),'PrettyPrint',true));
    fclose(fid);
else
    error('ORB5:TestAction','Use save or compare.');
end
disp(report);
end

function compare(a,b,path)
if isstruct(a)
    assert(isequal(fieldnames(a),fieldnames(b)),'Fields changed at %s.',path);
    fields = fieldnames(a);
    for k = 1:numel(fields), compare(a.(fields{k}),b.(fields{k}),[path '.' fields{k}]); end
elseif isnumeric(a)
    assert(isequal(size(a),size(b)),'Shape changed at %s.',path);
    assert(isequal(isnan(a),isnan(b)) && isequal(isinf(a),isinf(b)),'Nonfinite mask changed at %s.',path);
    finite = isfinite(a);
    scale = max(1,max(abs(double(a(finite)))));
    if isempty(scale), scale=1; end
    tolerance=1e-12;
    if isa(a,'single') || isa(b,'single'), tolerance=5e-6; end
    assert(all(abs(double(a(finite))-double(b(finite)))<=tolerance*scale),'Values changed at %s.',path);
else
    assert(isequaln(a,b),'Metadata changed at %s.',path);
end
end

function finish(visible)
close all force;
set(groot,'DefaultFigureVisible',visible);
end
