function run_orb5_tests
% Synthetic I/O + physics-regression + GUI lifecycle tests; no real run needed.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
folder=tempname(fullfile(root,'tests')); mkdir(folder);
cleanup=onCleanup(@()cleanupFolder(folder,root));
fid=fopen(fullfile(folder,'input'),'w');
fprintf(fid,['&equil\n a_mid=0.6, r0_mid=1.7\n btor0=2.0d0 ! tesla\n lx=720\n beta=0.01\n/\n' ...
    '&fields\n sfmin=0.1\n sfmax=0.9\n nsel_filter=''mn''\n nfilt1=0\n nfilt2=1\n n_flux_tube=1\n nchi=32\n/\n']);
fclose(fid);
assert(get_orb5_val('equil','btor0',fullfile(folder,'input'))==2);
assert(strcmp(get_orb5_val('fields','nsel_filter',fullfile(folder,'input')),'mn'));
file=fullfile(folder,'orb5_res.h5'); s=linspace(.1,.9,9).'; t=[0;2;5;9];
put(file,'/data/var1d/generic/flux_bin_psi',s.^2);
put(file,'/data/var0d/generic/time',t);
n=(2-s)*ones(1,4); T=(1-.3*s)*[1 1.1 1.2 1.3];
for name={'deuterium','fast'}
    base=['/data/var1d/' name{1} '/'];
    diagnostic(file,[base 'f_av'],n,t);
    diagnostic(file,[base 'v_par2_av'],1.5*n.*T,t);
    diagnostic(file,[base 'v_perp2_av'],1.5*n.*T,t);
    diagnostic(file,[base 'v_par_av'],zeros(size(n)),t);
    diagnostic(file,[base 'efluxw_rad'],ones(size(n)),t);
    diagnostic(file,[base 'pfluxw_rad'],2*ones(size(n)),t);
    put(file,['/data/var0d/' name{1} '/efluxw_tot'],t+1);
end
% An arbitrary additional ion with only one diagnostic and a different clock.
diagnostic(file,'/data/var1d/helium/efluxw_rad',3*ones(9,3),[1;4;8]);
for name={'deuterium','electrons','fast'}
    base=['/equil/profiles/' name{1} '/'];
    put(file,[base 's_prof'],s); put(file,[base 'rho_prof'],s);
    put(file,[base 'n_pic'],2-s); put(file,[base 't_pic'],1-s/3);
end
put(file,'/equil/profiles/generic/q',1+s.^2);
put(file,'/equil/profiles/generic/sgrid_eq',s);
config=orb5_config(struct('path',folder,'ad_hoc',0,'time_units','raw','visible','off'));
config.mass_ratios.fast=1; % Explicit test assumption, never inferred for an EP.
data=orb5_plot_library.load_flux_data(config);
assert(isequal(sort(data.species_names),sort({'deuterium','fast','helium'})));
assert(~isfield(data.species,'electrons'));
assert(isempty(data.f_ave));
assert(isfield(data.species.fast.diagnostics,'T'));
assert(~isfield(data.species.helium.diagnostics,'T')); % No guessed moments.
assert(max(abs(data.species.fast.diagnostics.T.values(:)-T(:)))<1e-10);
assert(isequal(data.species.helium.diagnostics.efluxw_rad.t,[1;4;8]));
eq=read_equil_data(file); assert(isfield(eq.species,'electrons'));
% Compound spectral storage uses the same dimension reversal as ORB5.
putPotential(file,9,12);
pot=orb5_plot_library.load_data(config);
assert(isequal(size(pot.p3d),[12 2 9 6]));
assert(~pot.is_em && isfield(pot,'EZ'));
assert(pot.t_norm(2)==1);
orb5_plot_library.plot_time_trace(pot,config,'phimax_lfs_Cs');
close all force;
theory=plot_dTi_theory(config,.5,'Verbose',false);
assert(all(theory.A(isfinite(theory.A))==0),'Default theory coefficient changed.');
close all force;
for mode={'radial_potential','temporal_potential','plot_potsc'}
    control=orb5_plot_library.cmd_interface(pot,config,mode{1});
    message=findall(control,'Tag','orb5_slice_status').Text;
    assert(startsWith(message,'Selected'),message);
    control.UserData.select(2);
    message=findall(control,'Tag','orb5_slice_status').Text;
    assert(startsWith(message,'Selected'),message);
    close(control);
end
% No electron transport is inferred from electron equilibrium profiles.
fig=orb5_plot_library.cmd_interface(data,config,'radial_flux');
controller=fig.UserData; controller.select(2);
figures=findall(groot,'Type','figure'); count=numel(figures);
controller.select(3);
assert(numel(findall(groot,'Type','figure'))==count,'Slice updates leaked plot windows.');
value=findall(fig,'Tag','orb5_slice_value'); assert(value.Value==data.t_norm(3));
exportapp(fig,fullfile(root,'tests','slice_preview.png'));
controller.select(1);
assert(numel(findall(groot,'Type','figure'))<=count);
close(fig); assert(~isvalid(controller.timer),'Controller timer leaked.');
singleData=data; singleData.s=single(data.s); singleData.t_norm=single(data.t_norm);
for mode={'radial_flux','temporal_flux'}
    control=orb5_plot_library.cmd_interface(singleData,config,mode{1});
    assert(startsWith(findall(control,'Tag','orb5_slice_status').Text,'Selected'));
    field=findall(control,'Tag','orb5_slice_value'); assert(isa(field.Value,'double'));
    if strcmp(mode{1},'radial_flux'), axisValues=singleData.t_norm; else, axisValues=singleData.s; end
    field.Value=double(axisValues(1)); field.ValueChangedFcn(field,[]);
    assert(startsWith(findall(control,'Tag','orb5_slice_status').Text,'Selected'));
    units=findall(control,'Type','uidropdown'); units.Value='Index'; units.ValueChangedFcn(units,[]);
    field.Value=numel(axisValues); field.ValueChangedFcn(field,[]);
    assert(field.Value==numel(axisValues));
    close(control);
end
app=ORB5_plotter(config);
app.UserData.inspect();
app.UserData.openPlot('equilibrium');
assert(startsWith(findall(app,'Tag','orb5_main_status').Value{1},'Opened'));
app.UserData.openPlot('calculator');
assert(startsWith(findall(app,'Tag','orb5_main_status').Value{1},'Opened'));
exportapp(app,fullfile(root,'tests','gui_preview.png'));
close(app); close all force;
% Cover EM preparation and the optional GAM/VMD branch with a small fixture.
putPotential(file,9,12,'pot3d_apar');
electromagnetic=config; electromagnetic.GAM=true;
em=orb5_plot_library.load_data(electromagnetic);
assert(em.is_em && isfield(em,'hfs_EM_max'));
assert(isequal(size(em.E_GAM),[12 9]));
orb5_plot_library.plot_time_trace(em,electromagnetic,'phimax_lfs_wA');
orb5_plot_library.plot_time_trace(em,electromagnetic,'radial_st');
close all force;
% No primary ion: fast and arbitrary species remain independently usable.
fid=H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
H5L.delete(fid,'/data/var1d/deuterium','H5P_DEFAULT');
H5L.delete(fid,'/data/var0d/deuterium','H5P_DEFAULT'); H5F.close(fid);
withoutIon=orb5_plot_library.load_flux_data(config);
assert(~isfield(withoutIon.species,'deuterium') && isfield(withoutIon.species,'fast'));
assert(isempty(withoutIon.f_avI));
% Scalar and vector harmonic extraction must agree, including nt=1.
obj=reshape(1:18,[1 1 3 6]); harms=reshape(repmat(0:5,3,1),[1 3 6]);
multi=get_m_ntor(obj,harms,[0 4]);
assert(isequal(reshape(multi(1,:,:),1,3),get_m_ntor(obj,harms,0)));
% Existing but incompatible data must fail, not disappear as optional.
diagnostic(file,'/data/var1d/helium/jdote',ones(7,3),[1;4;8]);
failed=false;
try, orb5_plot_library.load_flux_data(config);
catch err, failed=strcmp(err.identifier,'ORB5:Shape'); end
assert(failed,'Malformed existing diagnostic was silently accepted.');
fprintf('PASS: parser, discovered species, missing diagnostics, independent clocks, moments, GUI controls, timer cleanup, invalid shapes.\n');
end
function putPotential(file,ns,nt,name)
if nargin<4, name='pot3d'; end
base=['/data/var3d/generic/' name];
put(file,[base '/mmin'],zeros(ns,2)); put(file,[base '/time'],(0:nt-1).');
fid=H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
typ=H5T.create('H5T_COMPOUND',16);
H5T.insert(typ,'real',0,'H5T_NATIVE_DOUBLE');
H5T.insert(typ,'imaginary',8,'H5T_NATIVE_DOUBLE');
space=H5S.create_simple(4,[nt 2 ns 6],[]);
dataset=H5D.create(fid,[base '/data'],typ,space,'H5P_DEFAULT');
raw=reshape(1:6*ns*2*nt,[6 ns 2 nt]);
H5D.write(dataset,typ,'H5S_ALL','H5S_ALL','H5P_DEFAULT',struct('real',raw,'imaginary',raw/10));
H5D.close(dataset); H5S.close(space); H5T.close(typ); H5F.close(fid);
end
function put(file,path,value)
h5create(file,path,size(value)); h5write(file,path,value);
end
function diagnostic(file,base,value,t)
put(file,[base '/data'],value); put(file,[base '/time'],t);
end
function cleanupFolder(folder,root)
close all force;
if startsWith(folder,[fullfile(root,'tests') filesep]) && isfolder(folder), rmdir(folder,'s'); end
end
