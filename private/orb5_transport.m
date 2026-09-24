function varargout = orb5_transport(action, varargin)
    %ORB5_TRANSPORT Species diagnostics, transport preparation and plotting.
    %   data = orb5_transport('load', config)
    %   orb5_transport('plot', data, config, plot_id, slice_index)
    %
    % Every diagnostic owns its time axis. values(radius,time) is the canonical
    % layout; equilibrium-only species never acquire invented kinetic moments.
    switch action
        case 'load'
            varargout{1} = load_transport_data(varargin{:});
        case 'read'
            varargout{1} = read_flux_data(varargin{:});
        case 'plot'
            plot_transport_module(varargin{:});
        case 'convergence'
            plot_convergence(varargin{:});
        case 'theory'
            [varargout{1:nargout}] = plot_dTi_theory(varargin{:});
        otherwise
            error('ORB5:TransportAction','Unknown transport action: %s.',action);
    end
end

function data = load_transport_data(config)
    % Normalize each diagnostic's own clock before deriving moments.
    % Raw data and derived quantities share values(radius,time), t and s.
    config = orb5_config(config);
    data = read_flux_data(config.h5File);
    [factor,data.timelabel_str,data.plasma] = orb5_data('normalization', config);
    data.radiallabel_str = '$s$'; data.s = sqrt(double(data.psi_flux(:)));
    if config.ad_hoc && ~isempty(data.s)
        base = ['/equil/profiles/' config.profile_species '/'];
        data.s = interp1(h5read(config.h5File,[base 's_prof']), ...
            h5read(config.h5File,[base 'rho_prof']),data.s,'pchip');
        data.radiallabel_str = '$r/a$';
    end
    if any(~isfinite(data.s)) || any(diff(data.s)<=0)
        error('ORB5:Grid','Transport radial grid must be finite and strictly increasing.');
    end
    keys = fieldnames(data.species);
    % Collect clocks once to avoid repeatedly reallocating a growing vector.
    clocks = cell(1,numel(keys));
    for k = 1:numel(keys)
        key = keys{k}; sp = data.species.(key); ds = sp.diagnostics;
        quantities = fieldnames(ds);
        speciesClocks = cell(1,numel(quantities));
        for j = 1:numel(quantities)
            q = quantities{j}; d = ds.(q);
            if isempty(d.t), error('ORB5:Clock','Missing clock for %s/%s.',sp.name,q); end
            if any(~isfinite(d.t)) || any(diff(d.t)<=0)
                error('ORB5:Clock','Invalid clock for %s/%s.',sp.name,q);
            end
            d.t = double(d.t)/factor;
            if isempty(d.s)
                if numel(d.values)~=numel(d.t), error('ORB5:Shape','Invalid 0D diagnostic %s/%s.',sp.name,q); end
            else
                d.s = data.s;
                expected = [numel(d.s),numel(d.t)];
                if ~isequal(size(d.values),expected)
                    if isequal(size(d.values),fliplr(expected)), d.values = d.values.';
                    else, error('ORB5:Shape','Unexpected shape for %s/%s.',sp.name,q); end
                end
            end
            ds.(q) = d; speciesClocks{j} = d.t(:);
        end
        % Preserve the original moment expression, configurable mass and tau.
        required = {'f_av','v_par2_av','v_perp2_av','v_par_av'};
        if all(isfield(ds,required)) && isfield(config.mass_ratios,key)
            d = ds.f_av;
            compatible = true;
            for j = 2:numel(required)
                compatible = compatible && isequal(ds.(required{j}).t,d.t) ...
                    && isequal(size(ds.(required{j}).values),size(d.values));
            end
            if compatible && numel(d.s)>1
                n = d.values;
                pressure = config.mass_ratios.(key)/3 * ...
                    (ds.v_perp2_av.values + ds.v_par2_av.values - ds.v_par_av.values.^2);
                T = pressure./(n*config.tau); T(~isfinite(T)) = NaN;
                relativeT=(T-T(:,1))./T(:,1);
                win = min(11,size(T,1)); if mod(win,2)==0, win=win-1; end
                if strcmp(key,'deuterium') && win>=3
                    T = smoothdata(T,1,'sgolay',win);
                    n = smoothdata(n,1,'sgolay',win);
                end
                gradT = zeros(size(T)); gradLnT = gradT; gradLnn = gradT;
                for it = 1:size(T,2)
                    gradT(:,it) = gradient(T(:,it),d.s);
                    positiveT = T(:,it); positiveT(positiveT<=0)=NaN;
                    positiveN = n(:,it); positiveN(positiveN<=0)=NaN;
                    gradLnT(:,it) = -gradient(log(positiveT),d.s);
                    gradLnn(:,it) = -gradient(log(positiveN),d.s);
                end
                if strcmp(key,'deuterium'), w=5; else, w=11; end
                for it=1:size(T,2)
                    gradLnT(:,it) = smooth(gradLnT(:,it),min(w,size(T,1)));
                end
                if win>=3, gradLnn=smoothdata(gradLnn,1,'sgolay',win); end
                ratio = data.plasma.R0_cgs/data.plasma.a_cgs;
                ds.T = quantity(d,T);
                ds.RoLT = quantity(d,ratio*gradLnT);
                ds.RoLn = quantity(d,ratio*gradLnn);
                ds.dT = quantity(d,T-T(:,1));
                ds.dT_relative = quantity(d,relativeT);
                ds.dni_Ti = quantity(d,config.density_factor*(n-n(:,1)).*T(:,1));
                if isfield(ds,'efluxw_rad') && isequal(ds.efluxw_rad.t,d.t)
                    scale = (data.plasma.lx/2)^2;
                    if isfield(config.chi_scales,key) && ~isempty(config.chi_scales.(key))
                        scale = config.chi_scales.(key);
                    end
                    chi = -ds.efluxw_rad.values./(n.*gradT)*scale;
                    chi(~isfinite(chi)) = NaN;
                    ds.chi = quantity(d,chi);
                end
            end
        end
        clocks{k} = vertcat(speciesClocks{:});
        sp.diagnostics=ds; data.species.(key)=sp;
    end
    data.t_norm = unique(vertcat(clocks{:})); data.t = data.t_norm*factor;
    data.selection_t = data.t_norm;
    data.t_norm_0d = double(data.t0d)/factor;
    data = orb5_flux_aliases(data);
    data.t = data.t_norm*factor;
end
function d = quantity(template,values)
    d=template; d.values=values;
end


function data = orb5_flux_aliases(data)
    % Legacy fields for existing user scripts; new plots use species diagnostics.
    names = {'deuterium','electrons','fast'}; suffixes = {'I','e','E'};
    mapping = {'f_av','f_av';'f_av0','f_av';'df_av','df_av'; ...
        'pfluxw_rad','dpflux';'pfluxw_rad_neo','dpflux'; ...
        'efluxw_rad','deflux';'efluxw_rad_neo','deflux'; ...
        'efluxf0_rad','eflux';'efluxf0_rad_neo','eflux'; ...
        'v_par2_av','u2f';'v_perp2_av','vp2';'v_par_av','uf'; ...
        'efluxw_tot','deflux';'T','T';'chi','chi_perp';'RoLT','RoLT'};
    for k=1:3
        for j=1:size(mapping,1)
            q=mapping{j,1}; field=[mapping{j,2} suffixes{k}];
            if endsWith(q,'_neo'), field=[field '_neo']; end
            if strcmp(q,'f_av0'), field=[field '0']; end
            if strcmp(q,'efluxw_tot'), field=[field '_0d']; end
            if strcmp(q,'T') && k==1, field='Ti'; end
            if strcmp(q,'RoLT') && k==1, field='RoLTi'; end
            data.(field)=[];
            if isfield(data.species,names{k}) && isfield(data.species.(names{k}).diagnostics,q)
                data.(field)=data.species.(names{k}).diagnostics.(q).values;
            end
        end
    end
    if isfield(data.species,'deuterium')
        ds=data.species.deuterium.diagnostics;
        for pair={'dT','dTi';'dT_relative','d_Ti_Ti';'dni_Ti','dni_Ti';'RoLn','RoLn'}.'
            if isfield(ds,pair{1}), data.(pair{2})=ds.(pair{1}).values; end
        end
        if isfield(ds,'T')
            data.t_norm=ds.T.t;
        end
    end
end


function plot_transport_module(data,~,mode,index)
    % Species and quantities are independently optional; do not invent zeros.
    if nargin<4, index=1; end
    switch mode
        case {'radial_flux','temporal_flux'}
            quantities={'pfluxw_rad','pfluxw_rad_neo','efluxw_rad','efluxw_rad_neo','jdote'};
        case {'radial_transport','temporal_transport'}
            quantities={'chi','RoLn','RoLT','T'};
        case 'transport_st'
            quantities={'RoLT','chi','dT'};
        case 'transport_st_RoLT_component'
            quantities={'RoLT'};
        case 'eflux0d'
            quantities={'efluxw_tot'};
        case 'energy'
            quantities={'jdote_par','jdote_tot'};
        otherwise
            error('ORB5:Mode','Unknown transport mode: %s',mode);
    end
    keys=fieldnames(data.species); count=0;
    for j=1:numel(quantities)
        q=quantities{j}; ax=[];
        for k=1:numel(keys)
            sp=data.species.(keys{k}); ds=sp.diagnostics;
            if ~isfield(ds,q), continue; end
            d=ds.(q); if isempty(d.values), continue; end
            if startsWith(mode,'radial_')
                % Diagnostics may use different clocks. Select their nearest
                % stored time and report it; never extrapolate past coverage.
                target=data.t_norm(index);
                if target<min(d.t) || target>max(d.t), continue; end
                [~,ix]=min(abs(d.t-target)); x=d.s; y=d.values(:,ix);
                label=sprintf('%s (t=%.6g)',sp.name,d.t(ix));
                xlab=data.radiallabel_str;
            elseif startsWith(mode,'temporal_')
                target=data.s(index);
                [~,ix]=min(abs(d.s-target)); x=d.t; y=d.values(ix,:);
                label=sprintf('%s (r=%.6g)',sp.name,d.s(ix));
                xlab=data.timelabel_str;
            elseif startsWith(mode,'transport_st')
                f=figure('Name',[q ' : ' sp.name],'NumberTitle','off'); ax=axes(f);
                y=d.values;
                if strcmp(mode,'transport_st_RoLT_component'), y=y-y(:,1); end
                surface(ax,d.s,d.t,zeros(numel(d.t),numel(d.s)),y.', ...
                    'EdgeColor','none','FaceColor','interp'); view(ax,2); axis(ax,'tight');
                colorbar(ax); xlabel(ax,data.radiallabel_str,'Interpreter','latex');
                ylabel(ax,data.timelabel_str,'Interpreter','latex');
                title(ax,[strrep(mode,'_',' ') ' : ' sp.name ' / ' q],'Interpreter','none');
                count=count+1; continue;
            else
                x=d.t; y=d.values; label=sp.name; xlab=data.timelabel_str;
            end
            if isempty(ax)
                f=figure('Name',q,'NumberTitle','off'); ax=axes(f); hold(ax,'on');
            end
            plot(ax,x,y,'LineWidth',1.6,'DisplayName',label); count=count+1;
            xlabel(ax,xlab,'Interpreter','latex'); ylabel(ax,q,'Interpreter','none');
            legend(ax,'show','Interpreter','none'); grid(ax,'on');
            title(ax,strrep(mode,'_',' '));
        end
    end
    if count==0, error('ORB5:Unavailable','No available diagnostics for %s at this selection.',mode); end
end


function data = read_flux_data(h5File)
    %READ_FLUX_DATA Discover species and keep each diagnostic's own clock.
    paths = orb5_data('index', h5File);
    names = unique([orb5_data('species', paths,'/data/var1d'), ...
        orb5_data('species', paths,'/data/var0d')],'stable');
    data = struct('species',struct(),'species_names',{names});
    data.psi_flux = orb5_data('read_optional', h5File,paths,'/data/var1d/generic/flux_bin_psi');
    data.t0d = orb5_data('read_optional', h5File,paths,'/data/var0d/generic/time');
    rawRadius = sqrt(data.psi_flux(:));
    diagnostics = {'f_av','f_av0','df_av','pfluxw_rad','pfluxw_rad_neo', ...
        'efluxw_rad','efluxw_rad_neo','efluxf0_rad','efluxf0_rad_neo', ...
        'v_par2_av','v_par2_av0','v_perp2_av','v_perp2_av0', ...
        'v_par_av','v_par_av0','jdote'};
    for k = 1:numel(names)
        name = names{k}; key = matlab.lang.makeValidName(name);
        sp = struct('name',name,'diagnostics',struct());
        for j = 1:numel(diagnostics)
            q = diagnostics{j}; base = ['/data/var1d/' name '/' q];
            y = orb5_data('read_optional', h5File,paths,[base '/data']);
            if isempty(y), continue; end
            t = orb5_data('read_optional', h5File,paths,[base '/time']);
            if isempty(t)
                t = orb5_data('read_optional', h5File,paths,['/data/var1d/' name '/pfluxf0_rad/time']);
            end
            sp.diagnostics.(q) = struct('values',y,'t',t(:),'s',rawRadius);
        end
        for q = {'efluxw_tot','jdote_par','jdote_tot'}
            y = orb5_data('read_optional', h5File,paths,['/data/var0d/' name '/' q{1}]);
            if ~isempty(y)
                sp.diagnostics.(q{1}) = struct('values',y(:).','t',data.t0d(:),'s',[]);
            end
        end
        data.species.(key) = sp;
    end
end


function plot_convergence(data,config)
    %PLOT_CONVERGENCE Available species, each using its own diagnostic time grid.
    config=orb5_config(config); window=config.avg_window;
    validateattributes(window,{'numeric'},{'vector','numel',2,'finite'});
    if window(2)<window(1), error('ORB5:Window','Averaging window must increase.'); end
    keys=fieldnames(data.species); count=0;
    for q={'efluxw_rad','RoLT'}
        ax=[];
        for k=1:numel(keys)
            sp=data.species.(keys{k});
            if ~isfield(sp.diagnostics,q{1}), continue; end
            d=sp.diagnostics.(q{1}); mask=d.t>=window(1)&d.t<=window(2);
            if nnz(mask)<2, continue; end
            y=d.values(:,mask); avg=mean(y,2); spread=std(y,0,2);
            if isempty(ax), f=figure('Name',['Convergence_' q{1}]); ax=axes(f); hold(ax,'on'); end
            h=plot(ax,d.s,avg,'LineWidth',1.6,'DisplayName',sp.name);
            fill(ax,[d.s;flipud(d.s)],[avg+spread;flipud(avg-spread)],h.Color, ...
                'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off');
            xlabel(ax,data.radiallabel_str,'Interpreter','latex');
            ylabel(ax,q{1},'Interpreter','none'); legend(ax,'show','Interpreter','none'); grid(ax,'on');
            title(ax,sprintf('Mean and sample standard deviation, t = [%g, %g]',window));
            count=count+1;
        end
    end
    if count==0, error('ORB5:Window','No diagnostics have two samples in the requested averaging window.'); end
end


%% Experimental dTi comparison; preserve the implemented A and B definitions.
function result = plot_dTi_theory(config, r_target, varargin)
    % Experimental dTi comparison; preserve the implemented A and B definitions.
    config = orb5_config(config);
    C = config.theory_C; % Retain the original zero coefficient; scalar expands to the flux grid.
    if nargin < 2
        error('plot_dTi_theory:NotEnoughInputs', ...
            'Usage: plot_dTi_theory(config, r_target, ...).');
    end

    validateattributes(config, {'struct'}, {'nonempty'}, mfilename, 'config', 1);
    validateattributes(r_target, {'numeric'}, {'scalar', 'real', 'finite'}, ...
        mfilename, 'r_target', 2);

    parser = inputParser;
    parser.FunctionName = mfilename;
    addParameter(parser, 'Reload', false, ...
        @(x) islogical(x) || (isnumeric(x) && isscalar(x) && ismember(x, [0, 1])));
    addParameter(parser, 'C', C, ...
        @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));
    addParameter(parser, 'RadialMethod', 'interp', @(x) ischar(x) || isstring(x));
    addParameter(parser, 'TimeWindow', [], ...
        @(x) isempty(x) || (isnumeric(x) && numel(x) == 2 && diff(x) > 0));
    addParameter(parser, 'Verbose', true, ...
        @(x) islogical(x) || (isnumeric(x) && isscalar(x) && ismember(x, [0, 1])));
    addParameter(parser, 'FigurePosition', [0.35, 0.35, 0.34, 0.42], ...
        @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 4, 'finite'}));
    addParameter(parser, 'LineWidth', 2, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', '>', 0}));
    addParameter(parser, 'FontSize', 13, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', '>', 0}));
    addParameter(parser, 'LabelFontSize', 18, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', '>', 0}));
    addParameter(parser, 'LegendFontSize', 12, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', '>', 0}));
    addParameter(parser, 'AxesLineWidth', 2, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', '>', 0}));
    parse(parser, varargin{:});
    opts = parser.Results;

    opts.Reload = logical(opts.Reload);
    opts.Verbose = logical(opts.Verbose);
    opts.RadialMethod = validatestring(lower(char(opts.RadialMethod)), ...
        {'interp', 'nearest'}, mfilename, 'RadialMethod');
    config_local = local_normalize_config(config);
    persistent SESSION
    % This standalone comparison owns its prepared-data cache. Rebuild on a
    % changed session key or explicit Reload; changing only radius resamples it.

    need_reload = opts.Reload || isempty(SESSION) || ~isstruct(SESSION) || ...
        ~isfield(SESSION, 'key');

    if ~need_reload
        current_key = local_make_session_key(config_local);
        current_key.C = opts.C(:);          % reload if the C profile changed
        need_reload = ~isequaln(SESSION.key, current_key);
    end

    if need_reload
        if opts.Verbose
            fprintf(['plot_dTi_theory: loading ORB5 data ', ...
                '(first call, changed data/config/C, or Reload=true).\n']);
        end

        pot_data = orb5_plot_library.load_data(config_local);
        flux_data = orb5_plot_library.load_flux_data(config_local);

        session_key = local_make_session_key(config_local);
        session_key.C = opts.C(:);
        SESSION = local_build_session(pot_data, flux_data, config_local, ...
            session_key, opts);
        loaded_now = true;
    else
        if opts.Verbose
            fprintf(['plot_dTi_theory: reusing loaded session; ', ...
                'only the radial target is resampled.\n']);
        end
        loaded_now = false;
    end
    [t_out, y_A, y_B, y_AB, y_dTi, s_used, s_nearest] = ...
        local_sample(SESSION, r_target, opts);
    [fig1, ax1, h1] = local_plot_traces(SESSION, r_target, s_used, t_out, ...
        y_A, y_B, y_AB, y_dTi, opts);
    [fig2, ax2, cr, sr, ct, tt] = local_plot_correlation(SESSION, opts);
    result = struct( ...
        'r_target', r_target, ...
        'radial_method', opts.RadialMethod, ...
        'C', opts.C, ...
        'loaded_now', loaded_now, ...
        's_used', s_used, ...
        's_nearest', s_nearest, ...
        't', t_out(:), ...
        'A', y_A(:), ...
        'B', y_B(:), ...
        'AplusB', y_AB(:), ...
        'dTi', y_dTi(:), ...
        'corr_radial_s', sr(:), ...
        'corr_radial', cr(:), ...
        'corr_time', tt(:), ...
        'corr_temporal', ct(:), ...
        'figure_handle_traces', fig1, ...
        'axes_handle_traces', ax1, ...
        'figure_handle_corr', fig2, ...
        'axes_handle_corr', ax2, ...
        'line_handles_traces', h1);
end

function config_out = local_normalize_config(config_in)

    config_out = config_in;
    config_out = local_set_default(config_out, 'is_em', false);
    config_out = local_set_default(config_out, 'Z', 1);
    config_out = local_set_default(config_out, 'mu', 2);
    config_out = local_set_default(config_out, 'GAM', 0);
    config_out = local_set_default(config_out, 'smooth', 0);
    config_out = local_set_default(config_out, 'ad_hoc', 1);
    config_out.ad_hoc = 1;
end

function s = local_set_default(s, field_name, default_value)
    if ~isfield(s, field_name) || isempty(s.(field_name))
        s.(field_name) = default_value;
    end
end

function key = local_make_session_key(config)

    key = struct();
    key.pwd = config.path;
    key.config = config;

    key.files = struct( ...
        'p3dz_cache', local_file_stamp(fullfile(config.path, 'p3dz_cache.mat')), ...
        'flux_cache', local_file_stamp(fullfile(config.path, 'flux_cache.mat')), ...
        'orb5_res',   local_file_stamp(fullfile(config.path, 'orb5_res.h5')), ...
        'input',      local_file_stamp(fullfile(config.path, 'input')));
end

function stamp = local_file_stamp(file_name)
    d = dir(file_name);
    if isempty(d)
        stamp = struct('exists', false, 'datenum', 0, 'bytes', 0);
    else
        stamp = struct('exists', true, 'datenum', d(1).datenum, ...
            'bytes', d(1).bytes);
    end
end

function label = local_get_label(data, field_name, fallback)
    if isfield(data, field_name) && ~isempty(data.(field_name))
        label = data.(field_name);
    else
        label = fallback;
    end
end

function session = local_build_session(pot_data, flux_data, config_local, key, opts)

    session = struct();
    session.key = key;
    session.config = config_local;
    session.created = datetime('now');
    session.C = double(opts.C(:));      % C(r/a) profile as given (flux grid)
    session.pot_data = pot_data;
    session.flux_data = flux_data;

    session.timelabel_str = local_get_label(pot_data, 'timelabel_str', '$t$');
    session.radiallabel_str = local_get_label(pot_data, 'radiallabel_str', '$r/a$');
    if ~isfield(pot_data, 'EZZ') || isempty(pot_data.EZZ)
        error('plot_dTi_theory:MissingEZZ', ...
            ['pot_data.EZZ is missing or empty.  The requested quantity ', ...
            'requires the n = 0 zonal-flow field from the p3d data.']);
    end
    if ~isfield(pot_data, 't_norm') || isempty(pot_data.t_norm)
        error('plot_dTi_theory:MissingFieldTime', ...
            'pot_data.t_norm is missing or empty.');
    end
    [sE, EZZ] = local_prepare_time_major_axis(pot_data.s, pot_data.EZZ, 'EZZ');
    tE = local_prepare_axis_vector(pot_data.t_norm, [], 'field time');
    [tE, ~] = local_sort_axis(tE, 'field time');
    d2phi = -real(EZZ);                 % [nt x ns] on the field grid
    sE = double(sE);
    tE = double(tE);
    if ~isfield(flux_data, 't_norm') || isempty(flux_data.t_norm)
        error('plot_dTi_theory:MissingFluxTime', ...
            'flux_data.t_norm is missing or empty.');
    end
    if ~isfield(flux_data, 'dTi') || isempty(flux_data.dTi)
        error('plot_dTi_theory:MissingDTi', ...
            'flux_data.dTi is missing or empty (needs the flux data set).');
    end
    if ~isfield(flux_data, 'dni_Ti') || isempty(flux_data.dni_Ti)
        error('plot_dTi_theory:MissingDniTi', ...
            'flux_data.dni_Ti is missing or empty (needs the flux data set).');
    end
    if ~isfield(flux_data, 'Ti') || isempty(flux_data.Ti)
        error('plot_dTi_theory:MissingTiField', ...
            'flux_data.Ti is missing; needed for the initial temperature.');
    end

    sF_raw = flux_data.s(:).';
    [sF, order] = local_sort_axis(double(sF_raw), 'flux radial grid');
    tF_raw = flux_data.t_norm(:).';
    [tF, ~] = local_sort_axis(double(tF_raw), 'flux time');

    Ti0 = double(flux_data.Ti(:, 1));
    Ti0 = Ti0(order);
    dTi = double(flux_data.dTi);
    dTi = dTi(order, :);                % [ns x nt] on the flux grid
    B = -double(flux_data.dni_Ti);
    B = B(order, :);                    % [ns x nt] on the flux grid
    if isscalar(opts.C), opts.C=repmat(opts.C,numel(sF),1); end
    if numel(opts.C) ~= numel(sF)
        error('plot_dTi_theory:CSizeMismatch', ...
            ['The C profile has %d elements but the flux radial grid ', ...
            'has %d points. Set config.theory_C or the C option ', ...
            'to a scalar or a profile matching the flux grid.'], ...
            numel(opts.C), numel(sF));
    end
    C_field = interp1(sF, double(opts.C(:)), sE, 'linear');
    Ti0_field = interp1(sF, Ti0, sE, 'linear');
    % Preserve the implemented A definition: no division by Ti0 pending review.
    A = d2phi .* C_field;    % [nt x ns] on the field grid
    % Compare on the potential grid. Samples beyond flux coverage become NaN.
    F_B = griddedInterpolant({sF, tF}, B, 'linear', 'none');
    F_dTi = griddedInterpolant({sF, tF}, dTi, 'linear', 'none');
    B_field = F_B({sE, tE});          % [ns x nt] on the field grid
    dTi_field = F_dTi({sE, tE});      % [ns x nt] on the field grid

    session.sT = sE;                    % comparison radial grid = field grid
    session.tT = tE;                    % comparison time axis = field time
    session.tF = tF;                    % flux time range (correlation mask)
    session.Ti0 = Ti0;
    session.Ti0_field = Ti0_field;
    session.C_flux = double(opts.C(:)); % C(r/a) on the flux grid
    session.C_field = C_field;          % C(r/a) on the field grid
    session.A = A.';                    % [ns x nt]
    session.B = B_field;
    session.AplusB = A.' + B_field;
    session.dTi = dTi_field;
end

function [s_out, A_out] = local_prepare_time_major_axis(s_in, A_in, quantity_name)
    validateattributes(A_in, {'numeric'}, {'2d', 'nonempty'}, ...
        mfilename, quantity_name);
    s_out = local_prepare_axis_vector(s_in, size(A_in, 2), quantity_name);
    [s_out, order] = local_sort_axis(s_out, quantity_name);
    A_out = real(A_in);
    if ~isequal(order, 1:numel(s_out))
        A_out = A_out(:, order);
    end
end

function s_out = local_prepare_axis_vector(s_in, expected_length, quantity_name)
    validateattributes(s_in, {'numeric'}, {'vector', 'nonempty'}, ...
        mfilename, quantity_name);
    s_out = s_in(:).';
    if ~isempty(expected_length) && numel(s_out) ~= expected_length
        error('plot_dTi_theory:RadialGridSizeMismatch', ...
            '%s length (%d) does not match the data dimension (%d).', ...
            quantity_name, numel(s_out), expected_length);
    end
    if any(~isfinite(s_out))
        error('plot_dTi_theory:NonFiniteAxis', ...
            '%s contains NaN or Inf.', quantity_name);
    end
end

function [s_out, order] = local_sort_axis(s_in, quantity_name)
    if issorted(s_in)
        s_out = s_in;
        order = 1:numel(s_in);
    else
        [s_out, order] = sort(s_in);
    end

    if any(diff(s_out) <= 0)
        error('plot_dTi_theory:DuplicateAxis', ...
            ['%s is not strictly increasing after sorting; ', ...
            'cannot interpolate uniquely.'], quantity_name);
    end
end

function local_assert_inside_grid(r, s, quantity_name)
    if r < s(1) || r > s(end)
        error('plot_dTi_theory:RadiusOutsideGrid', ...
            ['Requested r/a = %.6g is outside the %s radial range ', ...
            '[%.6g, %.6g].'], r, quantity_name, s(1), s(end));
    end
end

function [t_out, y_A, y_B, y_AB, y_dTi, s_used, s_nearest] = ...
        local_sample(session, r_target, opts)
    local_assert_inside_grid(r_target, session.sT, 'field radial grid');

    [~, i_nearest] = min(abs(session.sT - r_target));
    s_nearest = session.sT(i_nearest);

    switch opts.RadialMethod
        case 'nearest'
            s_used = s_nearest;
            i = i_nearest;
            y_A = session.A(i, :);
            y_B = session.B(i, :);
            y_AB = session.AplusB(i, :);
            y_dTi = session.dTi(i, :);
        case 'interp'
            s_used = r_target;
            y_A = interp1(session.sT, session.A, r_target, 'linear');
            y_B = interp1(session.sT, session.B, r_target, 'linear');
            y_AB = interp1(session.sT, session.AplusB, r_target, 'linear');
            y_dTi = interp1(session.sT, session.dTi, r_target, 'linear');
    end

    t_out = session.tT(:);
    y_A = y_A(:);
    y_B = y_B(:);
    y_AB = y_AB(:);
    y_dTi = y_dTi(:);
end

function [fig, ax, h] = local_plot_traces(session, r_target, s_used, t, ...
        y_A, y_B, y_AB, y_dTi, opts)

    fig = figure('Name', sprintf('dTi theory vs numerical at r/a = %.4f', r_target), ...
        'Units', 'normalized', 'Position', opts.FigurePosition);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    h(1) = plot(ax, t, y_A, ...
        'LineWidth', opts.LineWidth, ...
        'Color', [0.00, 0.45, 0.74], ...
        'DisplayName', '$A = C\,\partial_r^2\phi_Z/T_{i0}$');
    h(2) = plot(ax, t, y_B, ...
        'LineWidth', opts.LineWidth, ...
        'Color', [0.47, 0.67, 0.19], ...
        'DisplayName', '$B = -\delta n_i\,T_{i0}$');
    h(3) = plot(ax, t, y_AB, ...
        'LineWidth', opts.LineWidth, ...
        'Color', [0.85, 0.10, 0.10], ...
        'DisplayName', '$A + B$');
    h(4) = plot(ax, t, y_dTi, ...
        'LineWidth', opts.LineWidth, ...
        'Color', 'k', 'LineStyle', '--', ...
        'DisplayName', '$\delta T_i/T_i$ (numerical)');

    yline(ax, 0, ':k');

    xlabel(ax, session.timelabel_str, ...
        'Interpreter', 'latex', 'FontSize', opts.LabelFontSize);
    ylabel(ax, '$\delta T_i/T_i$', ...
        'Interpreter', 'latex', 'FontSize', opts.LabelFontSize);
    c_used = interp1(session.sT, session.C_field, s_used, 'linear');
    title(ax, sprintf(['$r/a = %.4f$ (field grid), $C = %.3f$', ...
        ' | $A$: $r/a = %.4f$'], s_used, c_used, s_used), ...
        'Interpreter', 'latex', 'FontSize', opts.LabelFontSize);

    legend(ax, h, 'Interpreter', 'latex', ...
        'FontSize', opts.LegendFontSize, 'Box', 'off', 'Location', 'best');
    set(ax, 'LineWidth', opts.AxesLineWidth, 'FontSize', opts.FontSize);
    hold(ax, 'off');
end

function [mask, t_use] = local_time_mask(session, opts)
    t = session.tT(:);
    lo = t(1);
    hi = t(end);
    if ~isempty(session.tF)
        lo = max(lo, session.tF(1));
        hi = min(hi, session.tF(end));
    end
    if ~isempty(opts.TimeWindow)
        lo = max(lo, opts.TimeWindow(1));
        hi = min(hi, opts.TimeWindow(2));
    end
    mask = (t >= lo) & (t <= hi);
    t_use = t(mask);
end

function cc = local_corr(x, y)
    v = isfinite(x) & isfinite(y);
    if nnz(v) > 2 && std(x(v)) > 0 && std(y(v)) > 0
        c = corrcoef(x(v), y(v));
        cc = c(1, 2);
    else
        cc = NaN;
    end
end

function [fig, ax, cr, sr, ct, tt] = local_plot_correlation(session, opts)

    [mask, t_use] = local_time_mask(session, opts);
    if ~any(mask)
        warning('plot_dTi_theory:EmptyTimeWindow', ...
            ['The correlation time window is empty ', ...
            '(no overlap of field/flux time ranges).']);
        mask(:) = true;
        t_use = session.tT(:);
    end

    A_ = session.AplusB(:, mask);
    D_ = session.dTi(:, mask);
    s = session.sT(:);
    cr = NaN(numel(s), 1);
    for j = 1:numel(s)
        cr(j) = local_corr(A_(j, :), D_(j, :));
    end
    sr = s;
    ct = NaN(numel(t_use), 1);
    for k = 1:numel(t_use)
        ct(k) = local_corr(A_(:, k), D_(:, k));
    end
    tt = t_use;

    fig = figure('Name', 'Correlation of A+B with numerical dTi', ...
        'Units', 'normalized', 'Position', opts.FigurePosition);
    ax1 = subplot(2, 1, 1);
    hold(ax1, 'on');
    box(ax1, 'on');
    grid(ax1, 'on');
    plot(ax1, sr, cr, '.-', 'LineWidth', opts.LineWidth, ...
        'MarkerSize', 8, 'Color', [0.00, 0.45, 0.74]);
    yline(ax1, 0, ':k');
    xlabel(ax1, session.radiallabel_str, ...
        'Interpreter', 'latex', 'FontSize', opts.LabelFontSize);
    ylabel(ax1, '$\mathrm{corr}(A+B,\,\delta T_i/T_i)$', ...
        'Interpreter', 'latex', 'FontSize', opts.LabelFontSize);
    title(ax1, 'Radial correlation (time traces)', 'FontSize', opts.LabelFontSize);
    set(ax1, 'LineWidth', opts.AxesLineWidth, 'FontSize', opts.FontSize);
    ax2 = subplot(2, 1, 2);
    hold(ax2, 'on');
    box(ax2, 'on');
    grid(ax2, 'on');
    plot(ax2, tt, ct, '.-', 'LineWidth', opts.LineWidth, ...
        'MarkerSize', 8, 'Color', [0.85, 0.10, 0.10]);
    yline(ax2, 0, ':k');
    xlabel(ax2, session.timelabel_str, ...
        'Interpreter', 'latex', 'FontSize', opts.LabelFontSize);
    ylabel(ax2, '$\mathrm{corr}(A+B,\,\delta T_i/T_i)$', ...
        'Interpreter', 'latex', 'FontSize', opts.LabelFontSize);
    title(ax2, 'Temporal correlation (radial profiles)', 'FontSize', opts.LabelFontSize);
    set(ax2, 'LineWidth', opts.AxesLineWidth, 'FontSize', opts.FontSize);

    ax = [ax1, ax2];
end
