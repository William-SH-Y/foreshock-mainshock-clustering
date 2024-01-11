% Stackingdata_V6 subscript
% Plot foreshock and mainshock AE from Simulation 6 and 7

sim1 = 6;
sim2 = 7;

sz_select = 'rad';
f=figure;
[x,y,sz,c] = scatter_plot_AE(stack_FS_AE, sim1, sz_select,stack_MS_AE);
scatter(x, y, sz, c, 'filled');
xlim([-0.04,0.04]);
ylim([-0.02,0.02]);

caxis([-5e5, 0]);
cbar = colorbar;
set(get(cbar,'Title'), 'String','Time [steps]');
colormap("default")
xlabel('Along dip [m]');
ylabel('Along strike [m]');
title('Foreshock AE from Simulation 6');
set(gca,'FontName', 'Times');
set(gca,'FontSize',12);
exportgraphics(f,'fs_sim6_v2.jpeg','Resolution',300);
%%
%hold on 
f = figure;
[x,y,sz,c] = scatter_plot_AE(stack_FS_AE, sim2, sz_select,stack_MS_AE);
scatter(x, y, sz, c, 'filled','d');
xlim([-0.04,0.04]);
ylim([-0.02,0.02]);
cbar = colorbar;
caxis([-5e5, 0]);
set(get(cbar,'Title'), 'String','Time [steps]');
colormap("default")
xlabel('Along dip [m]');
ylabel('Along strike [m]');
title('Foreshock AE from Simulation 7');
set(gca,'FontName', 'Times');
set(gca,'FontSize',12);
exportgraphics(f,'fs_sim7_v2.jpeg','Resolution',300);
%hold off
%%
f=figure;
[x,y,sz,c] = scatter_plot_AE(stack_MS_AE, sim1, sz_select,stack_MS_AE);
scatter(x, y, sz, c, 'filled');
xlim([-0.04,0.04]);
ylim([-0.02,0.02]);
cbar = colorbar;
caxis([0, 7000]);
set(get(cbar,'Title'), 'String','Time [steps]');
colormap(flip(colormap("autumn")));
xlabel('Along dip [m]');
ylabel('Along strike [m]');
title('Mainshock AE from Simulation 6');
set(gca,'FontName', 'Times');
set(gca,'FontSize',12);
exportgraphics(f,'ms_sim6_v2.jpeg','Resolution',300);

%
f=figure;
%hold on 
[x,y,sz,c] = scatter_plot_AE(stack_MS_AE, sim2, sz_select, ...
    stack_MS_AE);
scatter(x, y, sz, c, 'filled','d');
xlim([-0.04,0.04]);
ylim([-0.02,0.02]);
caxis([0, 7000]);
cbar = colorbar;
colormap(flip(colormap("autumn")));
set(get(cbar,'Title'), 'String','Time [steps]');
xlabel('Along dip [m]');
ylabel('Along strike [m]');
title('Mainshock AE from Simulation 7');
set(gca,'FontName', 'Times');
set(gca,'FontSize',12);
exportgraphics(f,'ms_sim7_v2.jpeg','Resolution',300);
%hold off


function [x,y,sz,c] = scatter_plot_AE(df, sim_id, sz_select,...
    ms_df)
    sel = df(:,end)==sim_id;%select the simulation
    x = df(sel,4);
    y = df(sel,5);
    z = df(sel,6);
    fault_angle = pi/6;
    ms_trans = [sin(fault_angle), 0, -cos(fault_angle); ...
                0, 1, 0; ...
                cos(fault_angle), 0, sin(fault_angle)];
    new_coord = ms_trans * [x'; y'; z'];
    x = new_coord(1,:)';%dip
    y = new_coord(2,:)';%strike

    if strcmp(sz_select,'rad') == 1 
        sz = df(sel,7)*5000;%stack_FS_norm(sel,3) - min(stack_FS_norm(sel,3)) + 50;
    elseif strcmp(sz_select,'mag') == 1 
        sz = (df(sel,8)-min(df(sel,8))+0.0001)*100;
    else
        sz = 50;
    end

    onset = ms_df(ms_df(:,end)==sim_id, 2);
    onset = onset(1);
    c = (df(sel,2) - onset);    
end