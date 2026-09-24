function varargout = orb5_potential(action, varargin)
    %ORB5_POTENTIAL Potential loading, plotting, spectral tools and theory.
    %   data = orb5_potential('load', config)
    %   orb5_potential('plot', plot_id, data, config, slice_index)
    %
    % Array convention: p3d(time, toroidal_mode, radius, poloidal_slot).
    % Mode indices refer to the compact set selected from the HDF5 file.
    % Public helper files retain their established names and delegate here.
    switch action
        case 'load'
            varargout{1} = load_potential_data(varargin{:});
        case 'read'
            varargout{1} = read_pot3d_data(varargin{:});
        case 'plot'
            render_potential(varargin{:});
        case 'harmonic_map'
            [varargout{1:nargout}] = get_mlist(varargin{:});
        case 'harmonics'
            varargout{1} = get_m_ntor(varargin{:});
        case 'hfs'
            varargout{1} = hfs(varargin{:});
        case 'cwt'
            [varargout{1:nargout}] = cwt_complex(varargin{:});
        case 'vmd'
            [varargout{1:nargout}] = VMD(varargin{:});
        otherwise
            error('ORB5:PotentialAction','Unknown potential action: %s.',action);
    end
end

function render_potential(mode, data, config, index)
    % One routing table keeps all potential implementations in this file.
    switch mode
        case 'phimax_lfs_Cs'
            plot_phimax(data,config,mode);
        case 'radial_st'
            plot_radial_st(data,config);
        case 'radial_potential'
            plot_radial_m(data,config,index);
        case 'temporal_potential'
            plot_temporal_potential(data,config,index);
        case 'plot_potsc'
            plot_potsc(data,config,index);
        case 'phi_s_omega'
            plot_phi_s_omega(data,config,index);
        case 'f_Cs'
            plot_frequency_max(data,config,mode);
        case 'cross_correlation'
            plot_cross_correlation(data,config,index);
        otherwise
            error('ORB5:Mode','Unknown potential plot: %s.',mode);
    end
end

%% Prepare compact spectral fields and derived zonal quantities.
function data = load_potential_data(config)
    % Prepare compact spectral fields and derived zonal quantities.
    fprintf('Loading Potential Data...\n');
    config = orb5_config(config);
    [modeIndices,modeNumbers,modeBounds,modeStride]=orb5_mode_selection(config);
    data = read_pot3d_data(config.h5File,modeIndices);
    available_em = ~isempty(data.p3dz_apar);
    if isempty(config.is_em), config.is_em = available_em; end
    if config.is_em && ~available_em
        error('ORB5:MissingApar','EM was requested but pot3d_apar is absent.');
    end
    data.is_em = config.is_em;
    p3d_real      = permute((data.p3dz.real), [4,3,2,1]);
    data.p3dz.real = [];                              % free immediately
    p3d_imag      = permute((data.p3dz.imaginary), [4,3,2,1]);
    data.p3dz.imaginary = [];
    data.p3d      = complex(p3d_real, p3d_imag);
    clear p3d_real p3d_imag;

    data = rmfield(data, 'p3dz');
    if config.is_em
        p3d_apar_real      = permute((data.p3dz_apar.real), [4,3,2,1]);
        data.p3dz_apar.real = [];                              % free immediately
        p3d_apar_imag      = permute((data.p3dz_apar.imaginary), [4,3,2,1]);
        data.p3dz_apar.imaginary = [];
        data.p3d_apar      = complex(p3d_apar_real, p3d_apar_imag);
        clear p3d_apar_real p3d_apar_imag;
    end
    data.p3m = data.p3m.';

    data.t = double(data.p3t(:));

    data.sz = [size(data.p3d,1), size(data.p3d,2), size(data.p3d,3), size(data.p3d,4)];
    data.nt = data.sz(1);
    if data.nt<2 || data.sz(3)<2, error('ORB5:Grid','Potential analysis needs at least two times and radii.'); end
    if numel(data.t)~=data.nt || any(~isfinite(data.t)) || any(diff(data.t)<=0), error('ORB5:Clock','Invalid potential time axis.'); end
    data.dt = data.t(2)-data.t(1);

    sf_min = get_orb5_val('fields', 'sfmin', config.inputFile);
    sf_max = get_orb5_val('fields', 'sfmax', config.inputFile);
    data.ns = data.sz(3);
    data.s = linspace(sf_min, sf_max, data.ns);
    if config.ad_hoc == 1
        s_prof = h5read(config.h5File, ['/equil/profiles/' config.profile_species '/s_prof']);
        rho_prof = h5read(config.h5File, ['/equil/profiles/' config.profile_species '/rho_prof']);
        s = data.s;
        data.s = interp1(s_prof, rho_prof, s.', 'pchip');
        data.radiallabel_str = '$r/a$';
    else
        data.radiallabel_str = '$s$';
    end

    energy_s=data.efield_s;
    if config.ad_hoc && ~isempty(energy_s), energy_s=interp1(s_prof,rho_prof,energy_s,'pchip'); end
    if isempty(energy_s), energy_s=data.s; end
    if ~isempty(data.efield_non_zonal) && ~isempty(data.efield_zonal)
        data.efield_no_zf = trapz(energy_s,data.efield_non_zonal,1);
        data.efield_zf = trapz(energy_s,data.efield_zonal,1);
        data.efield_no_zf = data.efield_no_zf.';
        data.efield_zf = data.efield_zf.';
        data.efield_zf = interp1(data.efield_time, data.efield_zf, data.t, 'pchip');
        data.efield_no_zf = interp1(data.efield_time, data.efield_no_zf, data.t, 'pchip');
    else
        data.efield_no_zf = []; data.efield_zf = [];
    end

    [data.mlist, data.mharms_all] = get_mlist(data.sz, data.p3m);
    data.radial_envelope = sqrt(sum(abs(data.p3d).^2, 4));

    data.lfs_ES = sum(data.p3d, 4);
    data.lfs_ES_max = max(abs(data.lfs_ES), [], 3);
    % Keep legacy metadata fields, but all array indices use compact modes.
    data.nfilt1 = modeBounds(1);
    data.nfilt2 = modeBounds(2);
    if ~isempty(modeStride), data.n_flux_tube = modeStride; end
    data.nlist=modeNumbers;
    data.nlist_nft=modeNumbers;
    data.valid_idx=1:numel(modeNumbers);
    data.valid_nz_idx=find(modeNumbers~=0);
    data.is_zonal=modeNumbers==0;
    if any(data.is_zonal)
        validateattributes(config.zonal_slot,{'numeric'},{'scalar','integer','>=',1,'<=',data.sz(4)});
        data.phi_Z = (squeeze(data.p3d(:, data.is_zonal, :, config.zonal_slot)));
        phi_Z = data.phi_Z;
        % Columns are radial coordinates; retain the original zonal-field sign.
        [EZ,~] = gradient(phi_Z,data.s,1);

        EZ = -1*EZ;

        EZ_smooth = orb5_data('smooth', EZ,3,41,2);
        [EZZ,~] = gradient(EZ_smooth,data.s,1);
        data.EZ = EZ;

        data.EZZ =  -1*orb5_data('smooth', EZZ,3,21,2);
        data.EZ_max = max(abs(data.EZ), [], 2);
        if config.GAM
            [data.E_GAM, data.E_ZFZF] = Separate(data,config);
        end
    end
    if config.is_em
        data = precompute_em_data(data);
    end
    data.p = ORB5_calculator(config.Z, config.mu, config.inputFile);
    data.t_norm_wA = data.t / data.p.wci_wA0;
    data.t_norm_Cs = data.t / data.p.wci_Cs;
    if config.is_em == 1
        data.t_norm = data.t_norm_wA;
        data.timelabel_str = '$\omega_{\rm A0}t$';
        data.freqlabel_str = '$\omega/\omega_{\rm A0}$';

    else
        data.t_norm = data.t_norm_Cs;
        data.timelabel_str = '$(C_{\rm s}/R)t$';
        data.freqlabel_str = '$\omega/(C_{\rm s}/R)$';
    end
    if strcmp(config.time_units,'raw')
        data.t_norm=data.t; data.timelabel_str='$t$ (ORB5)'; data.freqlabel_str='$\omega$ (ORB5)';
    end
    data.dt_norm = data.t_norm(2) - data.t_norm(1);

    fprintf('Data preparation complete. Time range: [%.2f, %.2f]\n', ...
        min(data.t), max(data.t));
end


%% Resolve selected toroidal mode numbers to stored HDF5 slots.
function [indices,numbers,bounds,stride]=orb5_mode_selection(config)
    % Resolve selected toroidal mode numbers to stored HDF5 slots.
    info=h5info(config.h5File,'/data/var3d/generic/pot3d/data');
    count=info.Dataspace.Size(3);
    filter=get_orb5_val('fields','nsel_filter',config.inputFile);
    stride=[];
    if strcmp(filter,'mn')
        lo=get_orb5_val('fields','nfilt1',config.inputFile);
        hi=get_orb5_val('fields','nfilt2',config.inputFile);
        stride=get_orb5_val('fields','n_flux_tube',config.inputFile);
        numbers=lo:stride:hi;
    else
        numbers=get_orb5_val('fields','filtvaln',config.inputFile);
        lo=min(numbers); hi=max(numbers);
    end
    bounds=[lo hi];
    numbers=numbers(:).';
    if count==hi-lo+1
        indices=numbers-lo+1;
    elseif count==numel(numbers)
        indices=1:count;
    else
        error('ORB5:Modes','Input toroidal mode selection does not match the HDF5 mode dimension.');
    end
end


%% Maximum amplitudes, phase-frequency estimates and field energy.
function plot_phimax(data, config, ~)
    % Maximum amplitudes, phase-frequency estimates and field energy.
    if config.is_em == 1
        figure_name = 'phimax_lfs_wA';
    else
        figure_name = 'phimax_lfs_Cs';
    end

    fig_pos = [0.35, 0.35, 0.3, 0.4];

    figure('Name', figure_name,'Units', 'normalized', ...
        'Position', fig_pos);
    hold on; box on; grid on;
    if isfield(data,'EZ')
        EZ_max = data.EZ_max;
        [t_growth,growth_rate_est,fit_line] = get_growth_rate(EZ_max,data.t_norm,config);
        plot(data.t_norm, EZ_max, 'k', 'LineWidth', 3, ...
            'DisplayName', '$e\delta\partial_s\phi_{Z,max}/T_e$');
        plot(t_growth, fit_line, 'k', 'LineWidth', 3, ...
            'DisplayName', ['$fitting: \gamma=$', num2str(growth_rate_est, '%.4f')]);
    end
    for i = data.valid_nz_idx
        ntor = data.nlist(i);
        lfs_ES_max  = data.lfs_ES_max(:,i);

        lfs_ES_mid = data.lfs_ES(:,i,round(data.ns/2));
        [t_growth,growth_rate_est,fit_line] = get_growth_rate(lfs_ES_max,data.t_norm,config);
        frequency = get_frequency(lfs_ES_mid,data.t_norm,config);
        p1 = plot(data.t_norm, lfs_ES_max, 'LineWidth', 3, ...
            'DisplayName', sprintf('$e\\delta\\phi_{%d,max}/T_e$ and $\\omega = %.6f$', ntor,frequency));
        plot(t_growth, fit_line, 'LineWidth', 3, ...
            'DisplayName', ['$fitting: \gamma=$', num2str(growth_rate_est, '%.6f')]);
        if config.is_em
            plot(data.t_norm, data.p.coeff_A_to_eTepsi * data.hfs_EM_max(:, i), ...
                '--', 'LineWidth', 2, 'Color', get(p1, 'Color'), ...
                'DisplayName', sprintf('$e\\delta\\psi_{%d,max}/T_e$', ntor));
        end
    end

    set(gca, 'YScale', 'log', 'LineWidth', 2, 'FontSize', 13);
    xlabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Amplitude', 'FontSize', 18);

    legend('show', 'Interpreter', 'latex', 'FontSize', 18, 'Box', 'off', 'Location', 'best');

    if isempty(data.efield_no_zf), return; end
    figure('Name', 'Field energy evolution','Units', 'normalized', ...
        'Position', fig_pos);
    hold on; box on; grid on;
    if isfield(data,'EZ')
        [t_growth,growth_rate_est,fit_line] = get_growth_rate(data.efield_zf,data.t_norm,config);
        plot(data.t_norm, data.efield_zf, 'Color','r', 'LineWidth', 2, ...
            'DisplayName', '$E_{\rm zonal}$');
        plot(t_growth, fit_line, 'Color', 'b', 'LineWidth', 3, ...
            'DisplayName', ['$fitting: \gamma=$', num2str(growth_rate_est, '%.4f')]);
    end

    [t_growth,growth_rate_est,fit_line] = get_growth_rate(data.efield_no_zf,data.t_norm,config);
    p1 = plot(data.t_norm, data.efield_no_zf, 'Color', 'blue', 'LineWidth', 2, ...
        'DisplayName', sprintf('$E_{nonzonal}$'));
    plot(t_growth, fit_line, 'LineWidth', 3,  ...
        'DisplayName', ['$fitting: \gamma=$', num2str(growth_rate_est, '%.6f')]);
    if config.is_em
        plot(data.t_norm, data.p.coeff_A_to_eTepsi * data.hfs_EM_max(:, i), ...
            '--', 'LineWidth', 2, 'Color', get(p1, 'Color'), ...
            'DisplayName', sprintf('$e\\delta\\psi_{%d,max}/T_e$', ntor));
    end
    set(gca, 'YScale', 'log', 'LineWidth', 2, 'FontSize', 13);
    xlabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Amplitude', 'FontSize', 18);
    legend('show', 'Interpreter', 'latex', 'FontSize', 18, 'Box', 'off', 'Location', 'best');

end


%% Raw and normalized radius-time maps.
function plot_radial_st(data, config)
    % Raw and normalized radius-time maps.
    fig_pos = [0.35, 0.35, 0.3, 0.4];

    for imode = data.valid_nz_idx
        ntor = data.nlist(imode);
        radial_envelope = squeeze(data.radial_envelope(:,imode,:));

        figure('Name', sprintf('radial_st_n%d', ntor), 'Units', 'normalized', ...
            'Position', fig_pos);
        norm_envelope = radial_envelope ./ max(radial_envelope, [], 2);
        p = 1;  % Display exponent; 1 preserves linear normalized amplitudes.
        display_data = norm_envelope .^ p;

        pcolor(data.s, data.t_norm, display_data);
        shading('interp');
        colormap(turbo(256));
        c = colorbar;
        c.Ticks = [0, 0.1^p, 0.3^p, 0.5^p, 0.7^p, 1.0^p];
        c.TickLabels = {'0', '0.1', '0.3', '0.5', '0.7', '1.0'};
        c.Label.String = '$A/A_{max}$';
        c.Label.Interpreter = 'latex';
        c.Label.FontSize = 14;

        xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        ylabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        title(sprintf('Radial Envelope Normalized, $n=%d$', ntor), 'Interpreter', 'latex');
        set(gca, 'LineWidth', 2, 'FontSize', 13);
        text(0.02, 0.98, sprintf('Display: (A/A_{max})^{%.1f}', p), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1], 'FontSize', 10);

        figure('Name', sprintf('radial_st_n%d', ntor), 'Units', 'normalized', ...
            'Position', fig_pos);

        pcolor(data.s, data.t_norm, radial_envelope);
        shading('interp');
        colormap(turbo(256));
        c = colorbar;
        c.Label.Interpreter = 'latex';
        c.Label.FontSize = 14;

        xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        ylabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        title(sprintf('Radial Envelope Raw, $n=%d$', ntor), 'Interpreter', 'latex');
        set(gca, 'LineWidth', 2, 'FontSize', 13);
    end
    if isfield(data,'EZ')
        figure('Name', 'Zonal_Er_spatiotemporal', 'Units', 'normalized', ...
            'Position', fig_pos);

        EZ_real = real(data.EZ);
        norm_EZ = EZ_real ./ max(abs(EZ_real), [], 2);

        p = 1;
        pcolor(data.s, data.t_norm, sign(norm_EZ) .* abs(norm_EZ) .^ p);
        shading('interp');
        colormap(turbo(256));
        c = colorbar;
        c.Ticks = [-1, -0.5^p, 0, 0.5^p, 1];
        c.TickLabels = {'-1', '-0.5', '0', '0.5', '1'};
        c.Label.String = '$E_r/E_{r,max}$';
        c.Label.Interpreter = 'latex';

        xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        ylabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        title('Zonal $E_r$ Evolution Normalized', 'Interpreter', 'latex');
        set(gca, 'LineWidth', 2, 'FontSize', 13);

        figure('Name', 'Zonal_Er_spatiotemporal', 'Units', 'normalized', ...
            'Position', fig_pos);

        EZ_real = real(data.EZ);
        pcolor(data.s, data.t_norm, EZ_real);
        shading('interp');

        colormap(turbo(256));

        c = colorbar;
        c.Label.Interpreter = 'latex';

        xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        ylabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        title('Zonal $E_r$ Evolution Raw', 'Interpreter', 'latex');
        set(gca, 'LineWidth', 2, 'FontSize', 13);

        if config.GAM
            figure('Name', 'Zonal_Er_GAM_spatiotemporal', 'Units', 'normalized', ...
                'Position', fig_pos);

            E_GAM_real = real(data.E_GAM);

            pcolor(data.s, data.t_norm, E_GAM_real);
            shading('interp');

            colormap(turbo(256));

            c = colorbar;
            c.Label.Interpreter = 'latex';

            xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
            ylabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
            title('GAM $E_r$ Evolution Raw', 'Interpreter', 'latex');
            set(gca, 'LineWidth', 2, 'FontSize', 13);

            figure('Name', 'Zonal_Er_ZFZF_spatiotemporal', 'Units', 'normalized', ...
                'Position', fig_pos);

            E_ZFZF_real = real(data.E_ZFZF);

            pcolor(data.s, data.t_norm, E_ZFZF_real);
            shading('interp');
            colormap(turbo(256));
            c = colorbar;
            c.Label.Interpreter = 'latex';

            xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
            ylabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
            title('ZFZF $E_r$ Evolution Raw', 'Interpreter', 'latex');
            set(gca, 'LineWidth', 2, 'FontSize', 13);

            figure('Name', 'Zonal_Er_ZFZF_spatiotemporal_normalized', 'Units', 'normalized', ...
                'Position', fig_pos);

            E_ZFZF_real = real(data.E_ZFZF);
            norm_EZFZF = E_ZFZF_real ./ max(abs(E_ZFZF_real), [], 2);

            pcolor(data.s, data.t_norm, norm_EZFZF);
            shading('interp');
            colormap(turbo(256));
            c = colorbar;
            c.Label.Interpreter = 'latex';

            xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
            ylabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
            title('ZFZF $E_r$ Evolution Raw', 'Interpreter', 'latex');
            set(gca, 'LineWidth', 2, 'FontSize', 13);
        end
    end
end


%% Strongest poloidal harmonics at one time; optional theory overlay.
function plot_radial_m(data, config, ti)
    % Strongest poloidal harmonics at one time; optional theory overlay.
    fig_pos = [0.35, 0.35, 0.3, 0.4];

    valid_idx = data.valid_nz_idx;

    for imode = valid_idx
        ntor = data.nlist(imode);
        m_list = data.mlist{imode};
        m_num = length(m_list);

        p3d_slice = data.p3d(ti,imode,:,:);
        num_extract = round(m_num * 0.75);
        s_idx = round((m_num - num_extract)/2) + 1;

        figure('Name', sprintf('radial_m_n%d_t%.0f', ntor, data.t(ti)), ...
            'Units', 'normalized', 'Position', fig_pos);
        hold on; box on; grid on;
        all_amps = zeros(num_extract, length(data.s));
        max_amps = zeros(num_extract, 1);

        for idx = 1:num_extract
            mharm = m_list(s_idx + idx - 1);
            mamp = get_m_ntor(p3d_slice, data.mharms_all(imode,:,:), mharm);
            all_amps(idx, :) = (mamp);
            max_amps(idx) = max(abs(all_amps(idx, :)));
        end
        [sorted_amps, sort_idx] = sort(max_amps, 'descend');
        top_n = min(10, num_extract);
        top_idx = sort_idx(1:top_n);
        m_values = zeros(top_n, 1);
        for i = 1:top_n
            m_values(i) = m_list(s_idx + top_idx(i) - 1);
        end
        [~, m_sort_idx] = sort(m_values, 'descend');
        top_idx_sorted = top_idx(m_sort_idx);
        colors = hsv(top_n);
        for i = 1:top_n
            idx = top_idx_sorted(i);
            mharm = m_list(s_idx + idx - 1);
            line_color = colors(i, :);
            line_width = 1.2 + 2.3 * (max_amps(idx) / max(sorted_amps(1),eps));  % [1.2, 3.5]

            plot(data.s, abs(all_amps(idx, :)), ...
                'LineWidth', line_width, ...
                'Color', line_color, ...
                'DisplayName', sprintf('$m=%d$', mharm));
        end
        plot(data.s, squeeze(data.radial_envelope(ti,imode,:)), ...
            'LineWidth', 3, ...
            'Color', [0.8, 0.1, 0.1], ...
            'LineStyle', '--', ...
            'DisplayName', '$A(s)$');
        legend('show', 'Interpreter', 'latex', 'FontSize', 10, ...
            'Box', 'on', 'Location', 'best', 'NumColumns', 2, ...
            'Color', [1 1 1], 'EdgeColor', [0.5 0.5 0.5]);

        xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        ylabel('Amplitude', 'FontSize', 18);
        title(sprintf('$n=%d$, $t=%.2f$ (top 10 modes)', ntor, data.t_norm(ti)), ...
            'Interpreter', 'latex');
        set(gca, 'LineWidth', 2, 'FontSize', 13);
    end
    if isfield(data,'EZ')
        figure('Name', sprintf('radial_Zonal_Er_t%.0f', data.t_norm(ti)),'Units', 'normalized', ...
            'Position', fig_pos);
        if config.theory_overlay, [EZ_S,EZ_B] = get_theory(data,ti,config); end
        hold on;
        plot(data.s, real(smooth(data.EZ(ti,:),5)), ...
            'LineWidth', 1, 'DisplayName', '$E_{Z}$', 'Color', 'black','LineStyle','--');
        if config.theory_overlay && ~isempty(data.valid_nz_idx)
            plot(data.s,EZ_S, ...
                'LineWidth', 1, 'DisplayName', '$E_{Z,S}$', 'Color', 'black');
            plot(data.s,EZ_B, ...
                'LineWidth', 1, 'DisplayName', '$E_{Z,B}$', 'Color', 'red','LineStyle','--');
            plot(data.s,EZ_B.' + EZ_S, ...
                'LineWidth', 1, 'DisplayName', '$E_{Z,total}$', 'Color', 'red','LineStyle','--');
        end
        hold off;
        xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        ylabel('Amplitude','Interpreter', 'latex', 'FontSize', 18);
        legend('show', 'Interpreter', 'latex', 'FontSize', 12, 'Box', 'off');
        set(gca, 'LineWidth', 2, 'Box', 'on', 'FontSize', 13);
        grid on;
        title(sprintf('Zonal Flow, $t=%.2f$', data.t_norm(ti)), 'Interpreter', 'latex');
    end
end


%% Mode envelopes and the zonal field at one radius.
function plot_temporal_potential(data,~, si)
    % Mode envelopes and the zonal field at one radius.

    fig_pos = [0.35, 0.35, 0.3, 0.4];

    s = data.s;
    figure('Name', sprintf('temporal_n~=0_at_s=%.3f', s(si)), ...
        'Units', 'normalized', 'Position', fig_pos);
    hold on; box on; grid on;
    for imode = data.valid_nz_idx
        ntor = data.nlist(imode);
        plot(data.t_norm, smooth(data.radial_envelope(:,imode,si),10), ...
            'LineWidth', 2, ...
            'Color', 'r', ...
            'DisplayName', sprintf('n=%.0f',ntor));
    end

    legend('show', 'Interpreter', 'latex', 'FontSize', 15, ...
        'Box', 'on', 'Location', 'best', ...
        'Color', [1 1 1], 'EdgeColor', [0.5 0.5 0.5]);
    xlabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Amplitude', 'FontSize', 18);
    title(sprintf('Temporal evolution of $\\delta\\phi_n$ at s = %.3f', s(si)), ...
        'Interpreter', 'latex');
    set(gca, 'LineWidth', 2, 'FontSize', 13);

    if isfield(data,'EZ')
        figure('Name', sprintf('temporal_n=0_at_s=%.3f', s(si)), ...
            'Units', 'normalized', 'Position', fig_pos);
        hold on; box on; grid on;
        plot(data.t_norm, smooth(real(data.EZ(:,si)),10), ...
            'LineWidth', 2, ...
            'Color', 'r', ...
            'DisplayName', sprintf('n=0'));
    end

    legend('show', 'Interpreter', 'latex', 'FontSize', 15, ...
        'Box', 'on', 'Location', 'best', ...
        'Color', [1 1 1], 'EdgeColor', [0.5 0.5 0.5]);
    xlabel(data.timelabel_str, 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$E_Z$', 'FontSize', 18,'Interpreter', 'latex');
    title(sprintf('Temporal evolution of $E_Z$ at s = %.3f', s(si)), ...
        'Interpreter', 'latex');
    set(gca, 'LineWidth', 2, 'FontSize', 13);

end


%% Reconstruct poloidal sections and intensity-weighted radial wavenumber.
function plot_potsc(data,config, ti)
    % Reconstruct poloidal sections and intensity-weighted radial wavenumber.

    fig_pos = [0.35, 0.35, 0.3, 0.4];
    mharms = data.mharms_all;
    for imode = data.valid_nz_idx
        ntor = data.nlist(imode);
        p3dz_single = squeeze(data.p3d(ti,imode,:,:));
        nchi =  get_orb5_val('fields', 'nchi', config.inputFile);
        phi = 0;
        [potsc, chigrid] = create_potsc(p3dz_single, imode, nchi, phi, ntor, mharms);
        [S, C] = meshgrid(data.s, chigrid);  % both (nchi+1) x ns
        R0 = data.p.R0_cgs;
        a = data.p.a_cgs;
        R = R0 + a .* S .* cos(C);
        Z =       a .* S .* sin(C);
        Rn = (R - R0) ./ a;
        Zn =  Z       ./ a;

        figure('Name', sprintf('poloidal section ntor=%.0f', ntor), ...
            'Units', 'normalized', 'Position', fig_pos);
        hold on; box on; grid on;
        pcolor(Rn, Zn, real(potsc));
        shading interp;
        clim_val = max(abs(real(potsc(:))));
        clim([-max(clim_val,eps), max(clim_val,eps)]);
        colormap('jet');
        colorbar;
        xlabel('$(R-R_0)/a$','Interpreter','latex','FontSize',13,'FontName','Times New Roman');
        ylabel('$Z/a$','Interpreter','latex');
        title('Poloidal cross-section');

    end

    if isempty(data.valid_nz_idx), return; end
    [~, kr_avg] = get_kr(potsc, data.s, nchi);
    figure('Name', sprintf('Radial k_r at t=%.3f', data.t_norm(ti)), ...
        'Units', 'normalized', 'Position', fig_pos);
    hold on; box on; grid on;
    plot(data.s, smooth(kr_avg,30), ...
        'LineWidth', 2, ...
        'Color', 'r', ...
        'DisplayName', sprintf('n=%.0f',ntor));
    legend('show', 'Interpreter', 'latex', 'FontSize', 15, ...
        'Box', 'on', 'Location', 'best', ...
        'Color', [1 1 1], 'EdgeColor', [0.5 0.5 0.5]);
    xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Amplitude', 'FontSize', 18);
    title(sprintf('Radial dependence of  kr at t = %.3f', data.t_norm(ti)), ...
        'Interpreter', 'latex');
    set(gca, 'LineWidth', 2, 'FontSize', 13);

end


%% Radius-frequency map evaluated only at the requested time.
function plot_phi_s_omega(data,config,ti)
    % Radius-frequency map evaluated only at the requested time.
    t = data.t_norm;
    dt = t(2) - t(1);
    for imode = data.valid_nz_idx
        signal = squeeze(data.lfs_ES(:,imode,:));
        [C, f] = cwt_complex(signal, 1/dt, config.wname, ti);
        C_ti = reshape(C,numel(f),data.ns);
        C_ti(abs(C_ti) < max(abs(C_ti), [], 'all') * 0.6) = NaN;
        figure('Name', 'Time evolution','Units', 'normalized');
        hold on; box on; grid on;
        pcolor(data.s, f, abs(C_ti))
        colormap(turbo(256))
        shading interp;
        colorbar;
        xlabel(data.radiallabel_str, 'Interpreter', 'latex', 'FontSize', 18);
        ylabel(data.freqlabel_str, 'FontSize', 18,'Interpreter','latex');
        ylim([0 1])
        legend('show', 'Interpreter', 'latex', 'FontSize', 18, 'Box', 'off', 'Location', 'best');

        grid on;
    end

end


%% Peak wavelet frequency over time; original normalization.
function plot_frequency_max(data,config,~)
    % Peak wavelet frequency over time; original normalization.
    fig_pos = [0.35, 0.35, 0.3, 0.4];
    for imode = data.valid_idx
        ntor = data.nlist(imode);
        dt = data.t_norm(2) - data.t_norm(1);
        max_idx = round(data.ns/2);
        if ntor ==0
            signal = data.phi_Z(:,max_idx);
            [ac_component, ~] = VMD(signal, 1/dt);
            norm_signal = ac_component;
        else
            signal = squeeze(sum(data.p3d(:,imode,max_idx,:), 4));
            env = abs(hilbert(real(signal)));
            norm_signal = signal ./ max(env,eps);

        end

        [C, f] = cwt_complex(norm_signal, 1/dt, config.wname);
        f_parab = frequency_interpolation(C,f);

        figure('Name', 'Frequency analysis','Units', 'normalized', ...
            'Position', fig_pos);
        subplot(2,1,1);
        hold on;
        plot(data.t_norm, real(signal), 'Color', 'r','DisplayName','original (real)');
        plot(data.t_norm, real(norm_signal), 'Color', 'black','DisplayName','norm signal (real)');
        hold off;
        title(sprintf('Signal comparison n = %d',ntor));
        xlabel(data.timelabel_str); ylabel('Amplitude');
        legend('show', 'Interpreter', 'latex', 'FontSize', 13, 'Box', 'off');
        grid on;
        subplot(2,1,2);
        plot(data.t_norm, abs(f_parab), 'LineWidth', 2, 'Color', 'r');
        title('Evolution of Peak Frequency');
        xlabel(data.timelabel_str,'FontSize',13,Interpreter='latex'); ylabel(data.freqlabel_str,Interpreter='latex');
        grid on;
        legend('show', 'Interpreter', 'latex', 'FontSize', 12, 'Box', 'off');
        linkaxes(findall(gcf,'type','axes'), 'x');
        title(sprintf('Signal comparison n = %d',ntor));

    end
end


%% Intensity-curvature versus zonal-field delay at one radius.
function plot_cross_correlation(data,config,si)
    % Intensity-curvature versus zonal-field delay at one radius.
    if ~isfield(data,'EZ') || isempty(data.valid_nz_idx)
        error('ORB5:Unavailable','Cross-correlation needs zonal and non-zonal modes.');
    end
    t=data.t_norm(:); window=config.time_windows;
    idx=t>=window(1) & t<=window(2);
    if nnz(idx)<2
        error('ORB5:Window','Correlation window [%g,%g] needs at least two samples in [%g,%g].',window(1),window(2),t(1),t(end));
    end
    traces=zeros(numel(t),numel(data.valid_nz_idx));
    for k=1:numel(data.valid_nz_idx)
        imode=data.valid_nz_idx(k);
        % Filter all time columns together along radius. This is equivalent to
        % the original per-time loop, including the one-sided boundary gradients.
        intensity=double(abs(reshape(data.radial_envelope(:,imode,:),numel(t),[]).').^2);
        intensity=orb5_data('smooth',intensity,6,41,1);
        first=orb5_data('radial_gradient',intensity,data.s);
        second=orb5_data('radial_gradient',first,data.s);
        traces(:,k)=second(si,:).';
        a=traces(idx,k); b=real(data.EZ(idx,si));
        if ~any(abs(a)>0) || ~any(abs(b)>0)
            error('ORB5:Correlation','Cannot normalize correlation of an all-zero signal.');
        end
        [lag,~,lags,correlation]=calculate_nonlinear_delay(a,b,1/data.dt_norm);
        figure('Name',sprintf('Cross-correlation n=%d',data.nlist(imode)));
        plot(lags,correlation,'LineWidth',1.5); grid on;
        xlabel('Time lag (displayed time units)'); ylabel('Correlation coefficient');
        title(sprintf('n=%d, peak lag = %.6g',data.nlist(imode),lag));
        xline(0,'--r'); xline(lag,'--g');
    end
    figure('Name','Cross-correlation time traces');
    ax1=subplot(2,1,1); hold(ax1,'on');
    for k=1:numel(data.valid_nz_idx)
        plot(ax1,t,traces(:,k),'LineWidth',1.5, ...
            'DisplayName',sprintf('n=%d',data.nlist(data.valid_nz_idx(k))));
    end
    legend(ax1,'show'); grid(ax1,'on');
    ylabel(ax1,'Intensity curvature');
    title(ax1,sprintf('Selected radial coordinate = %.6g',data.s(si)));
    ax2=subplot(2,1,2);
    plot(ax2,t,real(data.EZ(:,si)),'LineWidth',1.5); grid(ax2,'on');
    ylabel(ax2,'$E_Z$','Interpreter','latex');
    xlabel(ax2,data.timelabel_str,'Interpreter','latex');
    linkaxes([ax1 ax2],'x');
end


%% EM high-field-side maxima; one column per selected toroidal mode.
function data = precompute_em_data(data)
    % EM high-field-side maxima; one column per selected toroidal mode.

    data.hfs_EM_max = zeros(size(data.p3d, 1), length(data.nlist), 'like', data.p3d);

    for i = data.valid_nz_idx
        temp_hfs = hfs(data.p3d_apar(:, i, :, :), data.p3m(i, :));
        data.hfs_EM_max(:, i) = max(abs(temp_hfs), [], 3);
    end
end


%% Original two-component VMD split along time at each radius.
function [E_GAM, E_ZFZF] = Separate(data,~)
    % Original two-component VMD split along time at each radius.
    E_GAM = zeros(data.nt,data.ns);
    E_ZFZF = zeros(data.nt,data.ns);
    for i = 1:data.ns
        [ac_component, dc_component] = VMD(data.EZ(:,i), 1/data.dt);
        E_GAM(:,i) = ac_component;
        E_ZFZF(:,i) = dc_component;
    end

end


%% Fit log amplitude in the original 0.1%-1% growth interval.
function [t_growth,growth_rate_est,fit_line] = get_growth_rate(signal,t_norm,~)
    % Fit log amplitude in the original 0.1%-1% growth interval.
    t = t_norm;
    idx = find(signal> 0.001 * max(signal) & signal < 0.01 * max(signal));
    idx = idx(isfinite(signal(idx)));
    if numel(idx)<2, t_growth=[]; growth_rate_est=NaN; fit_line=[]; return; end
    t_growth = t(idx)';
    s_growth = signal(idx)';
    poly_coeffs = polyfit(t_growth, log(s_growth), 1);
    growth_rate_est = poly_coeffs(1);
    intercept = poly_coeffs(2);
    fit_line = exp(intercept) * exp(growth_rate_est * t_growth);
end


%% Estimate phase slope on the second half of the time trace.
function frequency = get_frequency(signal,t_norm,~)
    % Estimate phase slope on the second half of the time trace.
    t = t_norm;
    idx = round(size(t,1)*1/2):1:size(t,1);
    t_growth = t(idx)';
    s_growth = signal(idx)';

    raw_phase = angle(s_growth);
    phase = unwrap(raw_phase);
    p = polyfit(t_growth, phase, 1);
    frequency = p(1);
end


%% Refine the wavelet ridge with a three-bin parabolic estimate.
function omega = frequency_interpolation(A,f)
    % Refine the wavelet ridge with a three-bin parabolic estimate.
    [~,max_i] = max(abs(A),[],1);
    sz = size(A);
    totalscal = sz(1);
    Nt = sz(2);
    omega = zeros(1, Nt);
    for col = 1:Nt
        idx = max_i(col);
        if idx == 1 || idx == totalscal
            omega(col) = f(idx);
            continue;
        end
        y = abs(A(idx-1:idx+1, col));
        x = f(idx-1:idx+1);
        dx = x(2) - x(1);
        num = y(1) - y(3); % Python 代码里的 y[0]-y[2]
        den = 2 * (y(1) - 2*y(2) + y(3));

        if den ~= 0
            offset = (num / den) * dx;
        else
            offset = 0;
        end
        omega(col) = x(2) + offset;
    end

end


%% Fourier reconstruction on a poloidal-angle by radius grid.
function [potsc, chigrid] = create_potsc(p3dz_single, imode, nchi, phi,ntor,mharms)
    % Fourier reconstruction on a poloidal-angle by radius grid.

    if nargin < 6 || isempty(phi),   phi   = 0.0;    end
    chigrid = 2*pi * (0:nchi)' / nchi;
    ns = size(p3dz_single, 1);
    potsc = zeros(nchi+1, ns, 'like', complex(0));
    amp = squeeze(p3dz_single(:, :));   % ns x mwidth
    m   = squeeze(mharms(imode, :, :));         % ns x mwidth
    phase = exp(1i * (ntor*phi + chigrid .* permute(m, [3 1 2])));
    potsc = potsc + sum(permute(amp, [3 1 2]) .* phase, 3);
end


%% Local and intensity-weighted radial wavenumber; mask weak amplitudes.
function [kr_local, kr_avg] = get_kr(potsc, s_grid, nchi)
    % Local and intensity-weighted radial wavenumber; mask weak amplitudes.
    s_grid  = s_grid(:)';                          % [1 x ns]
    ns      = numel(s_grid);
    chigrid = linspace(0, 2*pi, nchi+1)';          % [(nchi+1) x 1]

    d_potsc_ds = zeros(nchi+1, ns);
    for ichi = 1:(nchi+1)
        d_potsc_ds(ichi,:) = gradient(potsc(ichi,:), s_grid);
    end
    abs2      = abs(potsc).^2;                     % intensity |delta_U|^2
    integrand = conj(potsc) .* (-1i .* d_potsc_ds);% delta_U* (-i d/ds delta_U)
    amp_thresh      = 0.1 * max(abs(potsc(:)));
    mask            = abs(potsc) > amp_thresh;

    kr_local        = real(integrand) ./ (abs2 + eps);
    kr_local(~mask) = NaN;

    weight        = abs2;
    weight(~mask) = 0;

    integrand_m        = integrand;
    integrand_m(~mask) = 0;
    num_kr  = trapz(chigrid, integrand_m, 1);      % [1 x ns] complex
    denom   = trapz(chigrid, weight,      1);      % [1 x ns] real

    kr_avg  = real(num_kr) ./ (denom + eps);       % [1 x ns] Re<kr>

end


%% Experimental theory coefficients remain configurable and unchanged.
function [EZ_S,EZ_B] = get_theory(data,ti,config)
    % Experimental theory coefficients remain configurable and unchanged.

    s = data.s;
    lx = config.theory.lx;

    tau = config.theory.tau;
    n=config.theory.n;
    kappan = config.theory.kappan;
    a = data.p.a_cgs;
    R0 = data.p.R0_cgs;
    q = polyval(config.theory.q_polynomial,s);
    epsilon = s*(a/R0);
    etai =config.theory.etai;
    chi_Z = config.theory.chi_coefficient*q.^2./sqrt(epsilon);
    ktrhos = (n*q.'./((s+eps).'*0.5*lx));
    coeff_S = 1*(config.theory.frequency_ratio*(1+etai)-1).*(R0/a)./(chi_Z);
    coeff_B = 1*lx^2./(kappan*8*tau*chi_Z);

    for imode = data.valid_nz_idx
        p3d_slice = data.p3d(1:ti,imode,:,:);
        m_list = data.mlist{imode};
        m_num = length(m_list);
        A = squeeze(data.lfs_ES(ti,imode,:));
        A = abs(A).^2;

        order = 2;
        framelen = config.theory.smoothing_window;
        num_extract = round(m_num * 1);
        s_idx = round((m_num - num_extract)/2) + 1;
        e_idx = s_idx + num_extract - 1;
        m_request = m_list(s_idx : e_idx);
        mharms_all = data.mharms_all(imode, :, :);
        all_amps = get_m_ntor(p3d_slice, mharms_all, m_request);

        slice = all_amps(:, 1:ti, :);
        % slice is harmonic x time x radius: the third output differentiates radius.
        [~, ~,A1ds] = gradient(conj(slice), 1, 1, data.s);
        [~, ~,A1d] = gradient(slice, 1, 1, data.s);

        RS = 1i * ktrhos .* squeeze(sum(slice .* A1ds - conj(slice) .* A1d, 1)); % time x radius
        RS = sgolayfilt(double(RS),order,framelen,[],2);
        [RS,~] = gradient(RS,data.s,1);
        EZ_S = trapz(data.t_norm(1:ti),RS,1).*coeff_S.';
        EZ_S = sgolayfilt(double(EZ_S),order,framelen,[],2);

        A = sgolayfilt(double(A),order,framelen,[],1);
        A1d = gradient(A,data.s);

        A2d = gradient(A1d,data.s);

        EZ_B = coeff_B.*(A+1*(1+etai)*2*tau*(4/lx^2).*chi_Z.*A2d);
        EZ_B = sgolayfilt(double(EZ_B),order,framelen,[],1);
    end

end


%% Normalized cross-correlation; lag uses displayed time units.
function [lag_time, max_corr, lags, correlation] = calculate_nonlinear_delay(A, B, fs)
    % Normalized cross-correlation; lag uses displayed time units.
    A_prime = A;
    B_prime = B;
    [correlation, lag_indices] = xcorr(A_prime, B_prime, 'coeff');
    lags = lag_indices / fs;
    [max_corr, idx] = max(abs(correlation));
    lag_time = lags(idx);

end


%% Public read_pot3d_data implementation; see entry-point help for array layouts.
function data = read_pot3d_data(h5File, modeIndices)
    % Public read_pot3d_data implementation; see entry-point help for array layouts.
    if nargin<2, modeIndices=[]; end
    data = struct();
    index = orb5_data('index', h5File);

    paths = struct(...
        'p3dz', '/data/var3d/generic/pot3d/data', ...
        'p3dz_apar', '/data/var3d/generic/pot3d_apar/data', ...
        'p3m', '/data/var3d/generic/pot3d/mmin', ...
        'p3t', '/data/var3d/generic/pot3d/time', ...
        'efield_zonal', '/data/var1d/generic/efield0_1D/data', ...
        'efield_non_zonal', '/data/var1d/generic/efield_nozf_1D/data', ...
        'efield_s', '/data/var1d/generic/efield0_1D/coord1', ...
        'efield_time', '/data/var1d/generic/efield0_1D/time');

    fields = fieldnames(paths);
    for i = 1:length(fields)
        field = fields{i};
        path=paths.(field);
        if ~isempty(modeIndices) && ismember(field,{'p3dz','p3dz_apar'}) && ismember(path,index)
            info=h5info(h5File,path); shape=info.Dataspace.Size;
            validateattributes(modeIndices,{'numeric'},{'vector','integer','positive','<=',shape(3)});
            pieces=cell(1,numel(modeIndices));
            for k=1:numel(modeIndices)
                count=shape; count(3)=1;
                pieces{k}=h5read(h5File,path,[1 1 modeIndices(k) 1],count);
            end
            result=struct(); members=fieldnames(pieces{1});
            for k=1:numel(members)
                arrays=cellfun(@(p)p.(members{k}),pieces,'UniformOutput',false);
                result.(members{k})=cat(3,arrays{:});
            end
            data.(field)=result;
        else
            data.(field) = orb5_data('read_optional', h5File,index,path);
            if strcmp(field,'p3m') && ~isempty(modeIndices)
                data.p3m=data.p3m(:,modeIndices);
            end
        end
    end
end


%% Public get_mlist implementation; see entry-point help for array layouts.
function [mlist,mharms_all] = get_mlist(sz,p3m)
    % sz = [time,mode,radius,slot], p3m(mode,radius) stores the first slot.
    % Preserve the ORB5 sign convention for the physical poloidal harmonic.
    mharms_all = zeros(sz(2),sz(3),sz(4));
    mwidth = sz(4);
    moffsets = meshgrid(0:mwidth-1, 0:sz(3)-1);
    mlist = cell(sz(2), 1); % Pre-allocate a cell column

    for imode = 1:1:sz(2)

        mharms = -(double(p3m(imode,:).') + moffsets);
        mlist{imode} = unique(mharms).';
        mharms_all(imode,:,:) = mharms;
    end
end


%% Public get_m_ntor implementation; see entry-point help for array layouts.
function mamp=get_m_ntor(obj,mharms,mharm)
    % obj(time,1,radius,slot), mharms(1,radius,slot).
    % Output is (requested_harmonic,time,radius), or (time,radius) for one m.
    nt=size(obj,1); ns=size(obj,3); width=size(obj,4);
    if numel(mharms)~=ns*width, error('ORB5:Harmonics','Harmonic metadata shape mismatch.'); end
    mharms=reshape(mharms,1,ns,width);
    values=reshape(obj(:,1,:,:),nt,ns,width);
    result=zeros(numel(mharm),nt,ns,'like',obj);
    for k=1:numel(mharm)
        mask=mharms==mharm(k);
        result(k,:,:)=reshape(sum(values.*mask,3),1,nt,ns);
    end
    if isscalar(mharm), mamp=reshape(result,nt,ns); else, mamp=result; end
end


%% Public hfs implementation; see entry-point help for array layouts.
function n_fft = hfs(p3dz,p3m)
    % Reconstruct chi=pi by alternating the signs of adjacent poloidal slots.
    % p3dz(time,1,radius,slot) and p3m(1,radius) describe one toroidal mode.
    sign_factor = double(-(2 * mod(p3m, 2) - 1));

    p3dz_double = double(p3dz);
    sum_even = sum(p3dz_double(:, :, :, 1:2:end), 4);

    sum_odd = sum(p3dz_double(:, :, :, 2:2:end), 4);
    sign_factor = reshape(sign_factor, [1, 1, size(sign_factor, 2), 1]);
    diff = sum_even - sum_odd;

    n_fft = sign_factor .* diff;

end


%% Public cwt_complex implementation; see entry-point help for array layouts.
function [cfs, f, f_parab] = cwt_complex(A, fs, wname, timeIndices)
    % A(time,channel); cfs(frequency,selected_time,channel).
    % f is angular frequency in fs units. Preserve the legacy scale grid.
    persistent kernelCache
    A = double(A);

    if isvector(A)
        A = A(:);
    end

    [Nt, Nx] = size(A);
    if nargin<4, timeIndices=1:Nt; end
    validateattributes(timeIndices,{'numeric'},{'vector','integer','positive','<=',Nt});

    totalscal = 1024;

    f_c = centfrq(wname);
    cparam = 2*f_c*totalscal;
    a = totalscal:-1:1;
    scal = (cparam./a).';
    f = 2*pi*scal2frq(scal,wname,1/fs);
    cfs = complex(zeros(length(scal), numel(timeIndices), Nx));

    if nargin>=4 && numel(timeIndices)<Nt
        % Cache only the wavelet kernel, never simulation data. Repeated slider
        % updates reuse it. Bound retained memory to 128 MiB and one (Nt,wname) pair.
        if ~isempty(kernelCache) && kernelCache.nt==Nt && strcmp(kernelCache.wavelet,wname)
            response=kernelCache.response;
        else
            kernelCache=[];
            impulse=zeros(2*Nt-1,1); impulse(Nt)=1;
            response=cwt(impulse,scal,wname);
            if 16*numel(response)<=128*1024^2
                kernelCache=struct('nt',Nt,'wavelet',wname,'response',response);
            end
        end
        for k=1:numel(timeIndices)
            weights=response(:,Nt+timeIndices(k)-(1:Nt));
            % The legacy toolbox CWT conjugates a complex column internally.
            % Its centered impulse response also preserves boundary padding.
            cfs(:,k,:)=reshape(weights*conj(A),length(scal),1,Nx);
        end
    else
        for ix = 1:Nx
            column = cwt(A(:,ix), scal, wname);
            cfs(:,:,ix) = column(:,timeIndices);
        end
    end
    if nargout>2
        [~,idx] = max(abs(cfs),[],1);
        f_parab = reshape(f(idx),numel(timeIndices),Nx);
    end

end


%% Public VMD implementation; see entry-point help for array layouts.
function [ac_component, dc_component] = VMD(signal, ~)
    % Public VMD implementation; see entry-point help for array layouts.
    % Preserve the original two-mode split and fixed penalty pending review.
    [imfs, res] = vmd(real(signal), 'NumIMFs', 2, 'PenaltyFactor', 2000);
    ac_component = imfs(:, 1)+res;
    dc_component = sum(imfs(:, 2:end), 2);

end
