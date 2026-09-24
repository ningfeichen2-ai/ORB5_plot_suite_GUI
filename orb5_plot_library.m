classdef orb5_plot_library < handle
%ORB5_PLOT_LIBRARY Compatible public facade; implementation lives in modules.
methods(Static)
    function data=load_data(config)
        data=orb5_potential('load',orb5_config(config));
    end
    function data=load_flux_data(config)
        data=orb5_transport('load',orb5_config(config));
    end
    function data=load_equil_data(config)
        config=orb5_config(config); data=read_equil_data(config.h5File);
    end
    function plot_time_trace(data,config,mode,varargin)
        config=orb5_config(config); entry=orb5_plot_library.find_plot(mode,config);
        orb5_interface('plot',entry,data,config,1);
    end
    function fig=cmd_interface(data,config,plot_type)
        config=orb5_config(config); entry=orb5_plot_library.find_plot(plot_type,config);
        fig=orb5_interface('slice',data,config,entry);
    end
    function plot_equil(config)
        config=orb5_config(config);
        orb5_equilibrium('plot',read_equil_data(config.h5File),config);
    end
    function plot_radial_flux(data,config,index)
        orb5_transport('plot',data,orb5_config(config),'radial_flux',index);
    end
    function plot_temporal_flux(data,config,index)
        orb5_transport('plot',data,orb5_config(config),'temporal_flux',index);
    end
    function plot_radial_transport(data,config,index)
        orb5_transport('plot',data,orb5_config(config),'radial_transport',index);
    end
    function plot_temporal_transport(data,config,index)
        orb5_transport('plot',data,orb5_config(config),'temporal_transport',index);
    end
    function entry=find_plot(id,config)
        if ismember(id,{'phimax_lfs_wA','f_wA'})
            id=strrep(id,'wA','Cs');
        end
        entries=orb5_plot_registry(config); index=find(strcmp({entries.id},id),1);
        if isempty(index)
            if strcmp(id,'hysterisis')
                error('ORB5:MissingImplementation','The original suite calls plot_hysterisis but contains no implementation. Add one through orb5_plot_registry.');
            end
            error('ORB5:Mode','Unknown plot: %s.',id);
        end
        entry=entries(index);
    end
end
end
