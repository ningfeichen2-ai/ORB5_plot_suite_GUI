function varargout = plot_convergence(varargin)
%PLOT_CONVERGENCE Plot window averages and spread of transport quantities.
% Implementation is grouped in private/orb5_transport.m.
% plot_convergence(data,config) uses prepared transport data and config.avg_window.
[varargout{1:nargout}] = orb5_transport('convergence',varargin{:});
end
