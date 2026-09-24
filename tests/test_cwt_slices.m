function test_cwt_slices
% Compare sparse-time evaluation with the toolbox's full transform.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
rng(27);
for nt=[31 32 101]
    A=randn(nt,3)+1i*randn(nt,3);
    for wavelet={'cgau3','cmor1-1'}
        [full,f]=cwt_complex(A,12,wavelet{1});
        times=[1 ceil(nt/2) nt];
        [slices,fs]=cwt_complex(A,12,wavelet{1},times);
        expected=full(:,times,:);
        relative=max(abs(slices(:)-expected(:)))/max(abs(expected(:)));
        assert(relative<1e-10,'Sparse-time CWT differs from full CWT: %g.',relative);
        assert(isequal(f,fs));
        % Reuse the kernel with different data and sampling frequency. Only
        % the frequency axis depends on fs; no simulation values are cached.
        changed = 2*A+0.3i;
        [fullChanged,fChanged] = cwt_complex(changed,24,wavelet{1});
        [cached,fCached] = cwt_complex(changed,24,wavelet{1},times);
        expectedChanged = fullChanged(:,times,:);
        assert(max(abs(cached(:)-expectedChanged(:)))<1e-10*max(abs(expectedChanged(:))));
        assert(isequal(fCached,fChanged) && isequal(fChanged,2*f));
    end
end
fprintf('PASS: CWT slices match full complex CWT, including both boundaries.\n');
end
