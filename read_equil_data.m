function varargout = read_equil_data(varargin)
%READ_EQUIL_DATA Read available equilibrium profiles for each species.
% Implementation is grouped in private/orb5_equilibrium.m.
% data = read_equil_data(h5File). Profiles remain on each species' own grid;
% missing optional quantities are empty arrays rather than invented values.
[varargout{1:nargout}] = orb5_equilibrium('read',varargin{:});
end
