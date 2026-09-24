function varargout = read_flux_data(varargin)
%READ_FLUX_DATA Read species diagnostics with their individual time axes.
% Implementation is grouped in private/orb5_transport.m.
% data = read_flux_data(h5File). Each species.<name>.diagnostics.<quantity>
% has values(radius,time), its raw ORB5 clock t, and radial coordinate s.
% Use orb5_plot_library.load_flux_data(config) for normalization and moments.
[varargout{1:nargout}] = orb5_transport('read',varargin{:});
end
