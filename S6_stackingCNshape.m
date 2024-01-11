% Stackingdata_V6 subscript
% cluster based on shape of the CN segments

rat_cl_jl_cn_fs = zeros(max(stack_JL_FS(:,4)), 10);

for j = 2:length(t)
    %figure(f);f=f+1;
    for i = 1:max(stack_JL_FS(:,4))
        cond = (stack_JL_FS(:,4)==i) & ...
            stack_JL_FS(:,1)>=t(j) & ...
            stack_JL_FS(:,1)<=t(j-1);
%         subplot(3,4,i);
%         scatter(stack_JL_bcoord(cond,4),...
%             stack_JL_bcoord(cond,5));
%         xlim([-mv_W/2-0.005,mv_W/2+0.005]);
%         ylim([-mv_W/2-0.005,mv_W/2+0.005]);

        rat_cl_jl_cn_fs(i,j) = sum(cond)...
            /...
            sum(stack_JL_FS(:,1)>=t(j) & stack_JL_FS(:,1)<=t(j-1))...
            *100;
    end    
end

%
% plot the percentage of different cn segment shapes evolving with time
fig = figure(f);f=f+1;
for i = 1:height(rat_cl_jl_cn_fs)
    subplot(2,4,i);
    plot(-t(2:end), rat_cl_jl_cn_fs(i,2:end));
    xticks(-6e5:1e5:0);
    set(gca,'FontName', 'Times');
    %set(gca, 'XDir','reverse');
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Time [steps]');
ylabel('Proportion of data in each cluster [%]');
title('Normal force')
set(gca,'FontSize',12,'FontName', 'Times');

%exportgraphics(fig,'CN_cluster_proportion_evol.jpg','Resolution',300);

% cluster based on shape of the CS segments
rat_cl_jl_cs_fs = zeros(max(stack_JL_FS(:,4)), 10);
t = linspace(max(stack_JL_FS(:,1)),min(stack_JL_FS(:,1)),10);
for j = 2:length(t)
    %figure(f);f=f+1;
    for i = 1:max(stack_JL_FS(:,5))
        cond = (stack_JL_FS(:,5)==i) & ...
            stack_JL_FS(:,1)>=t(j) & ...
            stack_JL_FS(:,1)<=t(j-1);
%         subplot(3,4,i);
%         scatter(stack_JL_bcoord(cond,4),...
%             stack_JL_bcoord(cond,5));
%         xlim([-mv_W/2-0.005,mv_W/2+0.005]);
%         ylim([-mv_W/2-0.005,mv_W/2+0.005]);

        rat_cl_jl_cs_fs(i,j) = sum(cond)...
            /...
            sum(stack_JL_FS(:,1)>=t(j) & stack_JL_FS(:,1)<=t(j-1))...
            *100;
    end    
end
%
% plot the percentage of different cs segment shapes evolving with time
%fig = figure(f);f=f+1;
for i = 1:height(rat_cl_jl_cs_fs)
    subplot(2,4,i);
    plot(-t(2:end),rat_cl_jl_cs_fs(i,2:end));
    xticks(-6e5:1e5:0);
    set(gca,'FontName', 'Times');
    %set(gca, 'XDir','reverse');
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Time [steps]');
ylabel('Proportion of data in each cluster [%]');
title('Shear force')
set(gca,'Fontsize',12,'FontName', 'Times');
%exportgraphics(fig,'CS_cluster_proportion_evol.jpg','Resolution',300);

%%
% plot the percentage of different cn and cs segment shapes evolving with time
fig = figure(f);f=f+1;
legend_labels = {'normal force', 'shear force'};
for i = 1:height(rat_cl_jl_cn_fs)
    subplot(2,4,i);
    plot(-t(2:end), rat_cl_jl_cn_fs(i,2:end));
    hold on 
    plot(-t(2:end),rat_cl_jl_cs_fs(i,2:end));
    hold off
    xticks(-6e5:1e5:0);
    set(gca,'FontName', 'Times');
    %set(gca, 'XDir','reverse');
end
legend(legend_labels, 'Location', 'Best' ,'orientation','horizontal');
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Time [steps]');
ylabel('Proportion of data in each category [%]');
title('Temporal evolution of proportion');
set(gca,'FontSize',12,'FontName', 'Times');
exportgraphics(fig,'CNCS_cluster_proportion_evol.jpg','Resolution',300);
%%
% cluster based on shape of the CS segments
figure(f);f=f+1;
for i = 1:max(stack_JL_FS(:,5))
subplot(3,4,i);
scatter(stack_JL_bcoord(stack_JL_FS(:,5)==i,4),...
    stack_JL_bcoord(stack_JL_FS(:,5)==i,5));
xlim([-mv_W/2-0.005,mv_W/2+0.005]);
ylim([-mv_W/2-0.005,mv_W/2+0.005]);
end
