function entries = orb5_plot_registry(config)
%ORB5_PLOT_REGISTRY Add plots here, or append structs through config.extra_plots.
% Each entry has id, label, block, source, selector, callback(data,config,index).
if nargin==0, config=orb5_config(); end
entries = struct('id',{},'label',{},'block',{},'source',{},'selector',{},'callback',{});
add('phimax_lfs_Cs','Maximum amplitudes','Potential','potential','none',[]);
add('radial_st','Space-time envelopes','Potential','potential','none',[]);
add('radial_potential','Radial harmonics','Potential','potential','time',[]);
add('temporal_potential','Potential time traces','Potential','potential','radius',[]);
add('plot_potsc','Poloidal cross-section','Potential','potential','time',[]);
add('phi_s_omega','Radius-frequency spectrum','Potential','potential','time',[]);
add('f_Cs','Frequency evolution','Potential','potential','none',[]);
add('cross_correlation','Cross-correlation','Potential','potential','radius',[]);
add('radial_flux','Radial fluxes','Transport','transport','time',[]);
add('temporal_flux','Flux time traces','Transport','transport','radius',[]);
add('radial_transport','Radial transport','Transport','transport','time',[]);
add('temporal_transport','Transport time traces','Transport','transport','radius',[]);
add('transport_st','Transport space-time','Transport','transport','none',[]);
add('transport_st_RoLT_component','Change in R/LT','Transport','transport','none',[]);
add('eflux0d','Total energy flux','Transport','transport','none',[]);
add('energy','Power transfer J dot E','Transport','transport','none',[]);
add('convergence','Window averages and spread','Transport','transport','none',@(d,c,i)plot_convergence(d,c));
add('equilibrium','Equilibrium profiles','Equilibrium','equilibrium','none',[]);
add('calculator','Plasma parameters','ORB5 calculator','calculator','none',[]);
if ~isempty(config.extra_plots), entries=[entries config.extra_plots]; end
    function add(id,label,block,source,selector,callback)
        % Bind the plot id here; implementations remain local to the module.
        if strcmp(source,'potential')
            callback=@(d,c,i)orb5_potential('plot',id,d,c,i);
        end
        entries(end+1)=struct('id',id,'label',label,'block',block, ...
            'source',source,'selector',selector,'callback',callback);
    end
end
