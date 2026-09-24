function varargout = cwt_complex(varargin)
%CWT_COMPLEX Complex CWT, optionally evaluated only at selected times.
% Implementation is grouped in private/orb5_potential.m.
% [cfs,f,ridge] = cwt_complex(A,fs,wname[,timeIndices])
% A is time x channel; cfs is frequency x selected_time x channel.
% f and ridge are angular frequencies. Omitting timeIndices computes all times.
[varargout{1:nargout}] = orb5_potential('cwt',varargin{:});
end
