function varargout = VMD(varargin)
%VMD Variational mode decomposition; legacy interface retained.
% Implementation is grouped in private/orb5_potential.m.
% [ac,dc] = VMD(signal,fs). fs is retained for call compatibility;
% the existing toolbox decomposition operates on real(signal).
[varargout{1:nargout}] = orb5_potential('vmd',varargin{:});
end
