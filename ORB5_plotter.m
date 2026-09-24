function app=ORB5_plotter(config)
%ORB5_PLOTTER Open the GUI. ORB5_plotter(struct('path','run directory')).
% Keep only this version of the suite on the MATLAB path.
if nargin==0, config=struct(); end
config=orb5_config(config);
cache=struct(); entries=orb5_plot_registry(config);
app=uifigure('Name','ORB5 Plot Suite','Position',[120 120 1050 690],'Visible',config.visible);
main=uigridlayout(app,[4 1]); main.RowHeight={105,'1x',95,38};
settings=uigridlayout(main,[3 6]); settings.RowHeight={28 28 28};
uilabel(settings,'Text','Simulation folder');
folder=uieditfield(settings,'text','Value',config.path); folder.Layout.Column=[2 5];
uibutton(settings,'Text','Browse','ButtonPushedFcn',@browse);
uilabel(settings,'Text','Radial coordinate');
radial=uidropdown(settings,'Items',{'r/a','s'},'Value','r/a');
if ~config.ad_hoc, radial.Value='s'; end
uilabel(settings,'Text','Time coordinate');
time=uidropdown(settings,'Items',{'normalized','raw'},'Value',config.time_units);
uilabel(settings,'Text','Reference species');
profile=uieditfield(settings,'text','Value',config.profile_species);
uilabel(settings,'Text','Charge Z');
charge=uieditfield(settings,'numeric','Value',config.Z,'Limits',[eps Inf]);
uilabel(settings,'Text','Mass number');
mass=uieditfield(settings,'numeric','Value',config.mu,'Limits',[eps Inf]);
uibutton(settings,'Text','Reload / inspect','ButtonPushedFcn',@inspect);
uibutton(settings,'Text','Save config','ButtonPushedFcn',@saveConfig);
tabs=uitabgroup(main);
for block={'Potential','Transport','Equilibrium','ORB5 calculator'}
    tab=uitab(tabs,'Title',block{1}); grid=uigridlayout(tab,[2 1]); grid.RowHeight={'1x',36};
    subset=entries(strcmp({entries.block},block{1}));
    list=uilistbox(grid,'Items',{subset.label},'ItemsData',{subset.id});
    uibutton(grid,'Text','Open selected plot','ButtonPushedFcn',@(~,~)openPlot(list.Value));
end
status=uitextarea(main,'Editable','off','Tag','orb5_main_status','Value',{'Select a run folder containing input and orb5_res.h5.'; ...
    'Plots load on demand. Reload / inspect discards the in-memory session.'});
uilabel(main,'Text','Slice plots: slider + exact numeric entry + sample-index option + Play/Pause.');
app.UserData=struct('openPlot',@openPlot,'inspect',@inspect);
    function c=currentConfig()
        c=config; c.path=folder.Value; c.ad_hoc=strcmp(radial.Value,'r/a');
        c.time_units=time.Value; c.profile_species=profile.Value;
        c.Z=charge.Value; c.mu=mass.Value; c=orb5_config(c);
        if ~isfield(cache,'config') || ~isequaln(cache.config,c)
            cache=struct('config',c);
        end
    end
    function browse(~,~)
        result=uigetdir(folder.Value,'Choose ORB5 simulation folder');
        if isequal(result,0), return; end
        folder.Value=result; inspect();
    end
    function inspect(varargin)
        try
            c=currentConfig(); cache=struct('config',c);
            paths=orb5_data('index',c.h5File);
            eq=orb5_data('species',paths,'/equil/profiles');
            tr=unique([orb5_data('species',paths,'/data/var1d'),orb5_data('species',paths,'/data/var0d')]);
            status.Value={['Run: ' c.path]; ['Equilibrium species: ' strjoin(eq,', ')]; ...
                ['Transport species: ' strjoin(tr,', ')]; ...
                'Only available species diagnostics are plotted; files are read without using old MAT caches.'};
        catch err, status.Value={err.message}; end
    end
    function openPlot(id)
        try
            c=currentConfig(); entry=entries(strcmp({entries.id},id));
            status.Value={['Loading ' entry.label ' ...']}; drawnow;
            source=entry.source;
            if ~isfield(cache,source)
                switch source
                    case 'potential', cache.(source)=orb5_plot_library.load_data(c);
                    case 'transport', cache.(source)=orb5_plot_library.load_flux_data(c);
                    case 'equilibrium', cache.(source)=read_equil_data(c.h5File);
                    case 'calculator', cache.(source)=ORB5_calculator(c.Z,c.mu,c.inputFile);
                    otherwise, error('ORB5:Source','Unknown source: %s',source);
                end
            end
            data=cache.(source);
            if strcmp(entry.selector,'none')
                orb5_interface('redraw',@()orb5_interface('plot',entry,data,c,1),gobjects(0),c.visible);
            else
                orb5_interface('slice',data,c,entry);
            end
            status.Value={['Opened: ' entry.label]; 'Use Reload / inspect after simulation files change.'};
        catch err, status.Value={err.message}; end
    end
    function saveConfig(~,~)
        c=currentConfig(); [file,path]=uiputfile('*.mat','Save ORB5 configuration');
        if isequal(file,0), return; end
        config=c; save(fullfile(path,file),'config');
    end
end
