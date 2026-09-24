function test_correlation_vectorization
%TEST_CORRELATION_VECTORIZATION Compare plotted results with the original loop.
% Cover short smoothing windows, nonuniform radii and both radial boundaries.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
oldVisible = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
cleanup = onCleanup(@()finish(oldVisible));
rng(43);
config = orb5_config(struct('time_windows',[0 4]));
entry = orb5_plot_library.find_plot('cross_correlation',config);
for nr = [9 45]
    data.t_norm = (0:0.1:4).';
    data.dt_norm = 0.1;
    data.s = linspace(0.1,1,nr).^1.7;
    data.nlist = [0 4 8];
    data.valid_nz_idx = [2 3];
    data.timelabel_str = '$t$';
    data.EZ = single(randn(41,nr));
    data.radial_envelope = single(abs(randn(41,3,nr)));
    for si = [1 ceil(nr/2) nr]
        entry.callback(data,config,si);
        tracesFigure = findall(groot,'Type','figure','Name','Cross-correlation time traces');
        for mode = data.valid_nz_idx
            expected = zeros(41,1);
            for ti = 1:41
                intensity = double(abs(reshape(data.radial_envelope(ti,mode,:),[],1)).^2);
                window = min(41,nr-mod(nr+1,2));
                intensity = sgolayfilt(intensity,6,window);
                first = gradient(intensity,data.s);
                second = gradient(first,data.s);
                expected(ti) = second(si);
            end
            line = findall(tracesFigure,'Type','line', ...
                'DisplayName',sprintf('n=%d',data.nlist(mode)));
            assert(numel(line)==1);
            actual = line.YData(:);
            assert(max(abs(actual-expected))<1e-11*max(1,max(abs(expected))), ...
                'Vectorized curvature differs at radial index %d.',si);
            [expectedCorrelation,expectedLags] = xcorr(expected,real(data.EZ(:,si)),'coeff');
            correlationFigure = findall(groot,'Type','figure', ...
                'Name',sprintf('Cross-correlation n=%d',data.nlist(mode)));
            curve = findall(correlationFigure,'Type','line');
            assert(numel(curve)==1);
            assert(max(abs(curve.YData(:)-expectedCorrelation(:)))<1e-10);
            assert(max(abs(curve.XData(:)-expectedLags(:)*data.dt_norm))<1e-12);
        end
        close all force;
    end
end
fprintf('PASS: vectorized correlation matches the original loop at interior and boundary radii.\n');
end

function finish(visible)
close all force;
set(groot,'DefaultFigureVisible',visible);
end
