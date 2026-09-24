function varargout = read_pot3d_data(varargin)
%READ_POT3D_DATA Read selected potential modes; see private/orb5_potential.m.
% Implementation is grouped in private/orb5_potential.m.
% data = read_pot3d_data(h5File[,modeIndices]) reads raw compound HDF5 fields.
% modeIndices are one-based stored slots, not physical toroidal mode numbers.
% For prepared time x mode x radius x slot arrays, use load_data on the facade.
[varargout{1:nargout}] = orb5_potential('read',varargin{:});
end
