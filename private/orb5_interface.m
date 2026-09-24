function varargout = orb5_interface(action, varargin)
    %ORB5_INTERFACE Plot routing, slice controls and owned-window redraws.
    % Physics modules draw normal figures. This module manages their lifecycle;
    % it never stores another copy of the underlying simulation arrays.
    switch action
        case 'plot'
            orb5_dispatch(varargin{:});
        case 'slice'
            varargout{1} = orb5_slice_controls(varargin{:});
        case 'redraw'
            varargout{1} = orb5_redraw(varargin{:});
        otherwise
            error('ORB5:InterfaceAction','Unknown interface action: %s.',action);
    end
end

function orb5_dispatch(entry,data,config,index)
    config=orb5_config(config);
    if isfield(data,'is_em'), config.is_em=data.is_em; end
    if ismember(entry.id,{'f_Cs','phi_s_omega','cross_correlation'})
        dt=diff(data.t_norm);
        if isempty(dt) || max(abs(dt-dt(1)))>1e-6*abs(dt(1))
            error('ORB5:Sampling','This spectral/correlation plot requires uniformly spaced time samples.');
        end
    end
    if ~isempty(entry.callback)
        entry.callback(data,config,index); return;
    end
    if strcmp(entry.source,'transport')
        orb5_transport('plot', data,config,entry.id,index); return;
    elseif strcmp(entry.source,'equilibrium')
        orb5_equilibrium('plot', data,config); return;
    elseif strcmp(entry.source,'calculator')
        f=figure('Name','ORB5 calculator','NumberTitle','off');
        names=fieldnames(data); values=struct2cell(data);
        uitable(f,'Data',[names values],'ColumnName',{'Parameter (CGS unless named otherwise)','Value'}, ...
            'Units','normalized','Position',[0 0 1 1]); return;
    end
    error('ORB5:Callback','No callback for %s.',entry.id);
end


function fig=orb5_slice_controls(data,config,entry)
    if strcmp(entry.selector,'time')
        grid=double(data.t_norm(:)); axisLabel=data.timelabel_str;
    elseif strcmp(entry.selector,'radius')
        grid=double(data.s(:)); axisLabel=data.radiallabel_str;
    else
        error('ORB5:Selector','%s does not use a slice selector.',entry.id);
    end
    if isempty(grid) || any(~isfinite(grid)) || any(diff(grid)<=0)
        error('ORB5:Selector','The selection grid must be finite and strictly increasing.');
    end
    current=max(1,round(numel(grid)/2)); plots=gobjects(0); busy=false;
    fig=uifigure('Name',['ORB5 : ' entry.label],'Position',[200 180 650 240], ...
        'Visible',config.visible,'CloseRequestFcn',@closeController);
    layout=uigridlayout(fig,[4 5]); layout.RowHeight={28,44,30,'1x'};
    label=uilabel(layout,'Text',[entry.selector ' / ' axisLabel],'Interpreter','latex');
    label.Layout.Column=[1 2];
    units=uidropdown(layout,'Items',{'Coordinate','Index'},'ValueChangedFcn',@changeUnits);
    units.Layout.Column=3;
    value=uieditfield(layout,'numeric','Value',grid(current),'ValueChangedFcn',@typed);
    value.Layout.Column=[4 5]; value.Tag='orb5_slice_value';
    slider=uislider(layout,'Limits',[1 max(2,numel(grid))],'Value',current, ...
        'ValueChangedFcn',@slid); slider.Layout.Row=2; slider.Layout.Column=[1 5];
    slider.Tag='orb5_slice_slider';
    ticks=unique(round(linspace(1,numel(grid),min(6,numel(grid)))));
    slider.MajorTicks=ticks; slider.MinorTicks=[];
    slider.MajorTickLabels=arrayfun(@(i)sprintf('%.5g',grid(i)),ticks,'UniformOutput',false);
    if numel(grid)==1, slider.Enable='off'; end
    uibutton(layout,'Text','Previous','ButtonPushedFcn',@(~,~)select(current-1));
    uibutton(layout,'Text','Next','ButtonPushedFcn',@(~,~)select(current+1));
    play=uibutton(layout,'state','Text','Play','ValueChangedFcn',@playChanged);
    uilabel(layout,'Text','Step (samples)');
    step=uieditfield(layout,'numeric','Value',1,'Limits',[1 Inf],'RoundFractionalValues','on');
    status=uilabel(layout,'Text','','WordWrap','on','Tag','orb5_slice_status'); status.Layout.Row=4; status.Layout.Column=[1 5];
    timerObject=timer('ExecutionMode','fixedSpacing','Period',0.3,'BusyMode','drop','TimerFcn',@tick);
    fig.UserData=struct('timer',timerObject,'select',@select);
    select(current);
    function select(index)
        % Serialize slider/timer callbacks so redraws cannot overlap. Cleanup
        % releases the guard after success or an error in this selection.
        if busy || ~isvalid(fig), return; end
        busy=true; cleanup=onCleanup(@unlock);
        current=max(1,min(numel(grid),round(index))); slider.Value=current;
        if strcmp(units.Value,'Index'), value.Value=current; else, value.Value=grid(current); end
        try
            plots=orb5_redraw(@()orb5_dispatch(entry,data,config,current),plots,config.visible);
            status.Text=sprintf('Selected %s = %.9g; sample %d of %d.',entry.selector,grid(current),current,numel(grid));
        catch err
            status.Text=err.message; play.Value=false; stop(timerObject);
        end
    end
    function unlock, busy=false; end
    function typed(~,~)
        if strcmp(units.Value,'Index')
            if value.Value<1 || value.Value>numel(grid) || value.Value~=round(value.Value)
                status.Text='Enter an integer sample index inside the available range.'; return;
            end
            select(value.Value);
        else
            if value.Value<grid(1) || value.Value>grid(end)
                status.Text=sprintf('Enter a value between %.9g and %.9g.',grid(1),grid(end)); return;
            end
            [~,index]=min(abs(grid-value.Value)); select(index);
        end
    end
    function slid(~,event), select(event.Value); end
    function changeUnits(~,~)
        if strcmp(units.Value,'Index')
            value.Value=current;
            slider.MajorTickLabels=arrayfun(@num2str,ticks,'UniformOutput',false);
        else
            value.Value=grid(current);
            slider.MajorTickLabels=arrayfun(@(i)sprintf('%.5g',grid(i)),ticks,'UniformOutput',false);
        end
    end
    function playChanged(~,~)
        if play.Value, start(timerObject); else, stop(timerObject); end
    end
    function tick(~,~), select(mod(current-1+round(step.Value),numel(grid))+1); end
    function closeController(~,~)
        stop(timerObject); delete(timerObject);
        delete(plots(isgraphics(plots))); delete(fig);
    end
end


function output = orb5_redraw(callback,previous,visible)
    % Render offscreen, then replace contents of this controller's own figures.
    % A failed redraw leaves the last successful view intact.
    if nargin<2, previous=gobjects(0); end
    if nargin<3, visible='on'; end
    before=findall(groot,'Type','figure');
    old=get(groot,'DefaultFigureVisible');
    set(groot,'DefaultFigureVisible','off');
    cleanup=onCleanup(@()set(groot,'DefaultFigureVisible',old));
    try
        callback();
    catch err
        fresh=setdiff(findall(groot,'Type','figure'),before);
        delete(fresh); rethrow(err);
    end
    fresh=flipud(setdiff(findall(groot,'Type','figure'),before,'stable'));
    previous=previous(isgraphics(previous));
    output=gobjects(numel(fresh),1);
    for k=1:numel(fresh)
        if k<=numel(previous)
            target=previous(k); delete(target.Children);
            copyobj(fresh(k).Children,target);
            set(target,'Name',get(fresh(k),'Name'),'Visible',visible);
            delete(fresh(k)); output(k)=target;
        else
            output(k)=fresh(k); set(output(k),'Visible',visible);
        end
    end
    if numel(previous)>numel(fresh), delete(previous(numel(fresh)+1:end)); end
end
