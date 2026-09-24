function test_real_gui
% Exercise all main-window registered actions using the actual input files.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
config=orb5_config(struct('path',root,'visible','off'));
app=ORB5_plotter(config); cleanup=onCleanup(@()close('all','force'));
entries=orb5_plot_registry(config);
for entry=entries
    before=findall(groot,'Type','figure');
    app.UserData.openPlot(entry.id);
    message=findall(app,'Tag','orb5_main_status').Value;
    assert(startsWith(message{1},'Opened'),'%s: %s',entry.id,message{1});
    controls=findall(groot,'Tag','orb5_slice_status');
    for status=controls(:).'
        assert(startsWith(status.Text,'Selected'),'%s: %s',entry.id,status.Text);
    end
    fprintf('GUI PASS %s\n',entry.id);
    fresh=setdiff(findall(groot,'Type','figure'),before);
    for k=1:numel(fresh), if isgraphics(fresh(k)), close(fresh(k)); end; end
end
fprintf('PASS: all %d main GUI actions on the real run.\n',numel(entries));
end
