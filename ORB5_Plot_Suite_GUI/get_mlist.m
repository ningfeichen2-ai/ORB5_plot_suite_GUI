function varargout = get_mlist(varargin)
%GET_MLIST Map stored poloidal slots to physical harmonics.
% Implementation is grouped in private/orb5_potential.m.
% [mlist,mharms] = get_mlist(sz,p3m), sz = [time,mode,radius,slot].
% p3m(mode,radius) gives the first stored poloidal slot at each position.
[varargout{1:nargout}] = orb5_potential('harmonic_map',varargin{:});
end
