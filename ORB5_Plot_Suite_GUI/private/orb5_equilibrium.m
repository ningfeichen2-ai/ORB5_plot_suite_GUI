function varargout = orb5_equilibrium(action, varargin)
    %ORB5_EQUILIBRIUM Read and plot each species on its own equilibrium grid.
    % Stored gradient diagnostics keep their original derivative convention,
    % even when the horizontal coordinate is displayed as r/a.
    switch action
        case 'read'
            varargout{1} = read_equil_data(varargin{:});
        case 'plot'
            plot_equilibrium_module(varargin{:});
        otherwise
            error('ORB5:EquilibriumAction','Unknown equilibrium action: %s.',action);
    end
end

function data = read_equil_data(h5File)
    %READ_EQUIL_DATA Preserve each species' own radial coordinates.
    paths = orb5_data('index', h5File);
    names = orb5_data('species', paths,'/equil/profiles');
    data = struct('species',struct(),'species_names',{names});
    for k = 1:numel(names)
        sp = struct('name',names{k});
        for q = {'s_prof','rho_prof','n_pic','n_pic_mks','t_pic','t_pic_mks','gradn_pic','gradt_pic'}
            sp.(q{1}) = orb5_data('read_optional', h5File,paths,['/equil/profiles/' names{k} '/' q{1}]);
        end
        data.species.(matlab.lang.makeValidName(names{k})) = sp;
    end
    data.q = orb5_data('read_optional', h5File,paths,'/equil/profiles/generic/q');
    data.s = orb5_data('read_optional', h5File,paths,'/equil/profiles/generic/sgrid_eq');
end

function plot_equilibrium_module(data,config)
    keys=fieldnames(data.species); count=0;
    for quantity={'n_pic','t_pic','gradn_pic','gradt_pic'}
        q=quantity{1}; ax=[];
        for k=1:numel(keys)
            sp=data.species.(keys{k}); y=sp.(q);
            if config.ad_hoc, x=sp.rho_prof; label='$r/a$';
            else, x=sp.s_prof; label='$s$'; end
            if isempty(y) || isempty(x), continue; end
            if numel(y)~=numel(x), error('ORB5:Shape','Invalid equilibrium grid for %s/%s.',sp.name,q); end
            if startsWith(q,'grad'), y=-y; end
            if isempty(ax), f=figure('Name',q,'NumberTitle','off'); ax=axes(f); hold(ax,'on'); end
            plot(ax,x,y,'DisplayName',sp.name,'LineWidth',1.6); count=count+1;
            xlabel(ax,label,'Interpreter','latex');
            ylabel(ax,q,'Interpreter','none'); grid(ax,'on'); legend(ax,'show','Interpreter','none');
            if startsWith(q,'grad'), title(ax,['Minus stored ' q ' (original derivative coordinate)'],'Interpreter','none'); end
        end
    end
    if ~isempty(data.q) && ~isempty(data.s)
        x=data.s; label='$s$';
        if config.ad_hoc
            key=matlab.lang.makeValidName(config.profile_species);
            if ~isfield(data.species,key), error('ORB5:Profile','Missing reference species %s.',key); end
            sp=data.species.(key); x=interp1(sp.s_prof,sp.rho_prof,x,'pchip'); label='$r/a$';
        end
        f=figure('Name','Safety factor','NumberTitle','off'); ax=axes(f);
        plot(ax,x,data.q,'LineWidth',1.6); xlabel(ax,label,'Interpreter','latex'); ylabel(ax,'q'); grid(ax,'on');
        count=count+1;
    end
    if count==0, error('ORB5:Unavailable','No equilibrium profiles are available.'); end
end
