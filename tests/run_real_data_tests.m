function results=run_real_data_tests(sources)
% Exercise the real-file loaders, all registered plots, and slice callbacks.
% Does not write to the simulation input or HDF5 file.
if nargin<1, sources={'potential','transport','equilibrium','calculator'}; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
config=orb5_config(struct('path',root,'visible','off'));
output=fullfile(root,'tests','real_data'); if ~isfolder(output), mkdir(output); end
set(groot,'DefaultFigureVisible','off');
cleanup=onCleanup(@()close('all','force'));
entries=orb5_plot_registry(config); results=struct('id',{},'status',{},'seconds',{},'detail',{});
for source=sources
    source=source{1}; tic;
    try
        switch source
            case 'potential', data=orb5_plot_library.load_data(config);
            case 'transport', data=orb5_plot_library.load_flux_data(config);
            case 'equilibrium', data=orb5_plot_library.load_equil_data(config);
            case 'calculator', data=ORB5_calculator(config.Z,config.mu,config.inputFile);
        end
        record(['load_' source],'PASS',toc,'');
        if isfield(data,'t_norm')
            fprintf('TIME %s: [%g, %g]; %d samples\n',source,min(data.t_norm),max(data.t_norm),numel(data.t_norm));
        end
    catch err
        record(['load_' source],'FAIL',toc,getReport(err,'extended','hyperlinks','off'));
        continue;
    end
    subset=entries(strcmp({entries.source},source));
    for e=subset
        fprintf('START %s\n',e.id); start=tic;
        try
            if strcmp(e.selector,'none')
                if strcmp(source,'calculator')
                    app=ORB5_plotter(config); app.UserData.openPlot(e.id);
                    status=findall(app,'Tag','orb5_main_status');
                    assert(startsWith(status.Value{1},'Opened'),status.Value{1});
                    close(app);
                else
                    orb5_plot_library.plot_time_trace(data,config,e.id);
                end
            else
                controller=orb5_plot_library.cmd_interface(data,config,e.id);
                checkController(controller);
                if strcmp(e.selector,'time'), count=numel(data.t_norm); else, count=numel(data.s); end
                % A second, distinct slice catches window-update failures.
                controller.UserData.select(max(1,round(count/3)));
                checkController(controller);
            end
            drawnow;
            figures=findall(groot,'Type','figure');
            for k=1:numel(figures)
                if ~isempty(findall(figures(k),'Type','axes'))
                    exportgraphics(figures(k),fullfile(output,[e.id '.png']),'Resolution',100);
                    break;
                end
            end
            record(e.id,'PASS',toc(start),'');
        catch err
            record(e.id,'FAIL',toc(start),getReport(err,'extended','hyperlinks','off'));
        end
        close all force;
    end
    if strcmp(source,'potential')
        start=tic;
        try
            result=plot_dTi_theory(config,0.5,'Verbose',false);
            assert(~isempty(result.t));
            record('plot_dTi_theory','PASS',toc(start),'');
        catch err, record('plot_dTi_theory','FAIL',toc(start),getReport(err,'extended','hyperlinks','off')); end
        close all force;
        start=tic;
        try
            % The supplied run is electrostatic; test hfs reconstruction
            % numerically using its scalar spectral coefficients.
            mode=data.valid_nz_idx(1); block=data.p3d(1:3,mode,:,:);
            actual=hfs(block,data.p3m(mode,:));
            harmonics=reshape(data.mharms_all(mode,:,:),1,1,data.ns,[]);
            expected=sum(double(block).*((-1).^harmonics),4);
            assert(max(abs(actual(:)-expected(:)))<1e-10*max(1,max(abs(expected(:)))));
            record('hfs','PASS',toc(start),'Electrostatic coefficients used to verify reconstruction; no A-parallel dataset in this run.');
        catch err, record('hfs','FAIL',toc(start),getReport(err,'extended','hyperlinks','off')); end
        start=tic;
        try
            overlay=config; overlay.theory_overlay=true;
            entry=orb5_plot_library.find_plot('radial_potential',overlay);
            entry.callback(data,overlay,min(31,data.nt));
            record('get_theory_overlay','PASS',toc(start),'Runtime check at sample 31; physical model not revalidated.');
        catch err, record('get_theory_overlay','FAIL',toc(start),getReport(err,'extended','hyperlinks','off')); end
        close all force;
    end
    clear data;
end
    function record(id,status,seconds,detail)
        fprintf('%s %s (%.2fs)\n%s\n',status,id,seconds,detail);
        results(end+1)=struct('id',id,'status',status,'seconds',seconds,'detail',detail);
        fid=fopen(fullfile(output,['results_' strjoin(sources,'_') '.json']),'w');
        fprintf(fid,'%s',jsonencode(results,'PrettyPrint',true)); fclose(fid);
    end
end
function checkController(controller)
status=findall(controller,'Tag','orb5_slice_status');
assert(startsWith(status.Text,'Selected'),status.Text);
end
