%Stackingdata_V6 subscript
% Temporal evolution of particle speed


figure(f);f=f+1;
for j = 2:length(t)
    %subplot(3,4,j);
    %cond = stack_JL_FS(:,1)>=t(j) & ...
    %        stack_JL_FS(:,1)<=t(j-1); %within the time segment
    boxplot(stack_JL_FS(:,15)-stack_JL_FS(:,12), stack_JL_FS_timeseg);
end

figure(f);f=f+1;
for j = 2:length(t)
    boxplot(stack_JL_FS(:,14)-stack_JL_FS(:,11), stack_JL_FS_timeseg);
end

figure(f);f=f+1;
for j = 2:length(t)
    boxplot(...
        (stack_JL_FS(stack_JL_FS(:,10)>0,13)...
        -stack_JL_FS(stack_JL_FS(:,10)>0,10)) ...
        , ...
        stack_JL_FS_timeseg(stack_JL_FS(:,10)>0));
end

figure(f);f=f+1;
for j = 2:length(t)
    boxplot(...
        (stack_JL_FS(stack_JL_FS(:,10)<0,13)...
        -stack_JL_FS(stack_JL_FS(:,10)<0,10)    ) ...
        , ...
        stack_JL_FS_timeseg(stack_JL_FS(:,10)<0));
end
%histogram(stack_JL_FS(:,15)-stack_JL_FS(:,12));

figure(f);f=f+1;
for j = 2:length(t)
    boxplot(...
        stack_JL_FS(:,10),stack_JL_FS_timeseg...
        );
    ylabel('pre-foreshock particle velocity in x-dir');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    xlabel('from left to right: time windows approaching mainshocks')
end

figure(f);f=f+1;
for j = 2:length(t)
    boxplot(...
        stack_JL_FS(:,13),stack_JL_FS_timeseg...
        );
    ylabel('post-foreshock particle velocity in x-dir');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    xlabel('from left to right: time windows approaching mainshocks')
end
figure(f);f=f+1;
for j = 2:length(t)
    boxplot(...
        stack_JL_FS(:,11),stack_JL_FS_timeseg...
        );
    ylabel('pre-foreshock particle velocity in y-dir');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    xlabel('from left to right: time windows approaching mainshocks')
end
figure(f);f=f+1;
for j = 2:length(t)
    boxplot(...
        stack_JL_FS(:,14),stack_JL_FS_timeseg...
        );
    ylabel('post-foreshock particle velocity in y-dir');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    xlabel('from left to right: time windows approaching mainshocks')
end
figure(f);f=f+1;
for j = 2:length(t)
    boxplot(...
        stack_JL_FS(:,12),stack_JL_FS_timeseg...
        );
    ylabel('pre-foreshock particle velocity in z-dir');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    xlabel('from left to right: time windows approaching mainshocks')
end
figure(f);f=f+1;
for j = 2:length(t)
    boxplot(...
        stack_JL_FS(:,15),stack_JL_FS_timeseg...
        );
    ylabel('post-foreshock particle velocity in z-dir');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    xlabel('from left to right: time windows approaching mainshocks')
end