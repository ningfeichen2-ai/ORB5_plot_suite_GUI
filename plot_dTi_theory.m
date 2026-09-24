function varargout = plot_dTi_theory(varargin)
%PLOT_DTI_THEORY Plot the existing dTi theory comparison; use Reload=true after file changes.
% Implementation is grouped in private/orb5_transport.m.
% result = plot_dTi_theory(config,r_target,Name,Value,...).
% Options: Reload, C, RadialMethod ('interp'/'nearest'), TimeWindow, Verbose,
% FigurePosition, LineWidth, FontSize, LabelFontSize, LegendFontSize, AxesLineWidth.
% result includes sampled A/B/dTi traces, correlations and figure handles.
[varargout{1:nargout}] = orb5_transport('theory',varargin{:});
end
