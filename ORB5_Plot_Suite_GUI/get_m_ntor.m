function varargout = get_m_ntor(varargin)
%GET_M_NTOR Extract requested physical poloidal harmonics.
% Implementation is grouped in private/orb5_potential.m.
% amplitude = get_m_ntor(field,harmonicMap,requestedM).
% field is time x 1 x radius x slot, harmonicMap is 1 x radius x slot.
% Output is m x time x radius, or time x radius when requestedM is scalar.
[varargout{1:nargout}] = orb5_potential('harmonics',varargin{:});
end
