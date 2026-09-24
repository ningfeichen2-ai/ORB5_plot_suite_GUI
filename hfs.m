function varargout = hfs(varargin)
%HFS Reconstruct the potential on the high-field side.
% Implementation is grouped in private/orb5_potential.m.
% field = hfs(p3dz,p3m), with p3dz(time,1,radius,slot), p3m(1,radius).
[varargout{1:nargout}] = orb5_potential('hfs',varargin{:});
end
