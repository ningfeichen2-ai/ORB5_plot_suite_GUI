function varargout = orb5_data(action, varargin)
    %ORB5_DATA Shared HDF5 metadata, normalization and numerical utilities.
    % Missing optional datasets return []; existing unreadable datasets raise
    % errors, so file corruption is not mistaken for an absent particle species.
    switch action
        case 'index'
            varargout{1} = orb5_h5_index(varargin{:});
        case 'read_optional'
            varargout{1} = orb5_read_optional(varargin{:});
        case 'species'
            varargout{1} = orb5_species_names(varargin{:});
        case 'normalization'
            [varargout{1:nargout}] = orb5_normalization(varargin{:});
        case 'radial_gradient'
            varargout{1} = radial_gradient(varargin{:});
        case 'smooth'
            varargout{1} = orb5_sgolay(varargin{:});
        otherwise
            error('ORB5:DataAction','Unknown data action: %s.',action);
    end
end

function paths = orb5_h5_index(file)
    % Metadata index: absence is optional; corrupt existing data remains an error.
    paths = walk(h5info(file));
end
function paths = walk(group)
    % Allocate each subtree once rather than repeatedly growing the index.
    local = cell(numel(group.Datasets),1);
    for k = 1:numel(group.Datasets)
        local{k} = [regexprep(group.Name,'/$','') '/' group.Datasets(k).Name];
    end
    children = cell(numel(group.Groups),1);
    for k = 1:numel(group.Groups)
        children{k} = walk(group.Groups(k));
    end
    paths = [local; vertcat(children{:})];
end


function value = orb5_read_optional(file, paths, path)
    value = [];
    if ismember(path,paths), value = h5read(file,path); end
end


function names = orb5_species_names(paths,root)
    names = {};
    for k = 1:numel(paths)
        token = regexp(paths{k},['^' root '/([^/]+)/'],'tokens','once');
        if ~isempty(token) && ~strcmp(token{1},'generic')
            names{end+1} = token{1}; %#ok<AGROW>
        end
    end
    names = unique(names,'stable');
end


function [factor,label,p] = orb5_normalization(config)
    config = orb5_config(config);
    p = ORB5_calculator(config.Z,config.mu,config.inputFile);
    if strcmp(config.time_units,'raw')
        factor = 1; label = '$t$ (ORB5)';
    else
        em = config.is_em;
        if isempty(em)
            paths = orb5_h5_index(config.h5File);
            em = ismember('/data/var3d/generic/pot3d_apar/data',paths);
        end
        if em, factor = p.wci_wA0; label = '$\omega_{A0}t$';
        else, factor = p.wci_Cs; label = '$(C_s/R)t$'; end
        if ~isfinite(factor) || factor<=0
            error('ORB5:Normalization','Selected time normalization is not finite and positive.');
        end
    end
end


function y = orb5_sgolay(x,order,window,dim)
    % Reduce to an odd window that fits the selected dimension; very short
    % grids pass through unchanged rather than failing in sgolayfilt.
    n = size(x,dim); window = min(window,n-mod(n+1,2));
    if window <= order, y = x; return; end
    y = sgolayfilt(double(x),order,window,[],dim);
end

function result = radial_gradient(values, radius)
    %RADIAL_GRADIENT Differentiate each column along a shared radial grid.
    % Match MATLAB gradient's central differences and one-sided endpoints on
    % nonuniform grids. Keeping input arithmetic preserves its numeric precision.
    radius = radius(:);
    if size(values,1) ~= numel(radius) || numel(radius)<2
        error('ORB5:Grid','Radial differentiation needs at least two matching radii.');
    end
    result = zeros(size(values),'like',values);
    result(1,:) = (values(2,:)-values(1,:))/(radius(2)-radius(1));
    result(end,:) = (values(end,:)-values(end-1,:))/(radius(end)-radius(end-1));
    result(2:end-1,:) = (values(3:end,:)-values(1:end-2,:))./(radius(3:end)-radius(1:end-2));
end
