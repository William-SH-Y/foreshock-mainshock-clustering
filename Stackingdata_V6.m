% Update: Hausdorff distances calculated from Mahalanobis distance
cd 'D:\William Google Drive\Matlab Codes 2\Matlab Code'
clear
close all
clc
% importing
f = 1;
mv_W = 40e-3;
mv_H = 107e-3;

% farray = {'201','208','214','215','216','220','221','222','223','224','225', '226','227','228',...
%      '229','230','231','232','233','234','235', '236','237','238','239','240','241'};%,'242', '243', '244',...
%  %    '245','246','247','248','249','250','251','254','255'};%
farray = {'487','488','489','490','486','485','484','476','475','474','473','472',...
    '464','483','471','449','470','455','469','458','468','467','479','482','481','480','463'};%, '463-krat5'};
fgroup = num2str(400);
for i = 1:length(farray)
    load(strcat('S',farray{i},'jump_cat'),'jump_cat');
    jump_Cat(i) = jump_cat;
    %jump over jump_coef(slope of lin fit to)
    load(strcat('S',farray{i},'ratio_jump'),'ratio_jump');
    fnR = fieldnames(ratio_jump);
    jump_Rat(i) = ratio_jump;
    
    load(strcat('S',farray{i},'coef_jump'),'coef_jump');
    fnC = fieldnames(coef_jump);
    jump_Coef(i) = coef_jump;
    
    load(strcat('S',farray{i},'ladder_jump'),'ladder_jump');
    fn = fieldnames(ladder_jump);
    for k=1:numel(fn)
        if( isnumeric(ladder_jump.(fn{k})) )
            ladder_jump.(fn{k}) = ladder_jump.(fn{k})';
        end
    end
    jump_Ladder(i) = ladder_jump;
    
    
    load(strcat('S',farray{i},'sbs'),'sbs');
    Sbs{i} = sbs;
    
    load(strcat('S',farray{i},'bcoord'),'bcoord');
    Bcoord{i} = bcoord;
    Ccoord{i} = [(bcoord(:,1)+bcoord(:,4))/2 , (bcoord(:,2)+bcoord(:,5))/2, ...
        (bcoord(:,3)+bcoord(:,6))/2];
    
    load(strcat('S',farray{i},'seg_sjcrk'),'seg_sjcrk');
    seg_Sjcrk{i} = seg_sjcrk;
    
    load(strcat('S',farray{i},'sj_ss_ini'),'sj_stickslip_ini');
    ss_Ini{i} = sj_stickslip_ini;
    
    load(strcat('S',farray{i},'sj_kernel2000'),'ySix');
    sj_Kernel2000{i} = ySix;
    
    load(strcat('S',farray{i},'kernelpks'),'kernelpeaks');
    kn_Pks{i} = kernelpeaks;
    
    load(strcat('S',farray{i},'kernelpks_leverarm'),'kernelpks_leverarm');
    kn_leverarm{i} = kernelpks_leverarm;
    
    load(strcat('S',farray{i},'sj_spatdist'),'sj_spatdist');
    sj_spacdist{i} = sj_spatdist;
    
    load(strcat('S',farray{i},'pf_ss_sj'),'pf_ss_sj');
    sj_ss_polyfit{i} = pf_ss_sj;
    
    load(strcat('S',farray{i},'ss_stressdrop'),'ss_stressdrop');
    ss_sd{i} = ss_stressdrop; %stickslip stressdrop
    
    load(strcat('S',farray{i},'seg_ms_sj'),'seg_ms_sj');
    seg_AE{i} = seg_ms_sj; %stickslip stressdrop
    
%     load(strcat('S',farray{i},'broken_trigger'),'broken_trigger');
%     broken_Trigger{i} = broken_trigger; %stickslip stressdrop
%     load(strcat('S',farray{i},'broken_stresstrigger'),'broken_stresstrigger');
%     broken_ST{i} = broken_stresstrigger; %stickslip stressdrop
%     load(strcat('S',farray{i},'broken_inttriggerx'),'broken_inttriggerx');
%     broken_IXT{i} = broken_inttriggerx; %stickslip stressdrop
%     load(strcat('S',farray{i},'broken_inttriggery'),'broken_inttriggery');
%     broken_IYT{i} = broken_inttriggery; %stickslip stressdrop
%     load(strcat('S',farray{i},'broken_inttriggerz'),'broken_inttriggerz');
%     broken_IZT{i} = broken_inttriggerz; %stickslip stressdrop
    
    load(strcat('S',farray{i},'seg_cn_prejump'),'seg_cn_prejump');
    seg_CN_pre{i} = seg_cn_prejump; %stickslip stressdrop
    load(strcat('S',farray{i},'seg_cs_prejump'),'seg_cs_prejump');
    seg_CS_pre{i} = seg_cs_prejump; %stickslip stressdrop
    
    load(strcat('S',farray{i},'cyc_wrtss'),'cyc_wrtss');
    cyc_Wrtss{i} = cyc_wrtss; 
    load(strcat('S',farray{i},'bd_mean'),'bd_mean');
    bd_Mean{i} = bd_mean; 
    load(strcat('S',farray{i},'bd_size'),'bd_size');
    bd_Size(i,:) = bd_size; 
end

fault_angle = pi/6;
coord_trans = [sin(fault_angle), 0, -cos(fault_angle); ...
    0, 1, 0; ...
    cos(fault_angle), 0, sin(fault_angle)];

for i = 1:length(farray)
    temp = zeros(size(seg_Sjcrk{i},1), 1);
    for j = 1:size(seg_Sjcrk{i},1)
        temp(j) = size(seg_Sjcrk{i}{j},1);
    end
    window_stickslip(i) = find(temp == max(temp));
    clear temp
end


% Stacking foreshocks
L = 0;
for i = 1:length(farray)
    for j = 1:window_stickslip(i)-1
        L = L + length(seg_AE{i}{j}(:,8));
    end
end
W = length(seg_AE{1}{1}(1,:));

stack_FS_AE = zeros(L, W+1);
stack_FS_time = zeros(L, 1);
k = 1;
for i = 1:length(farray)
    for j = 1:window_stickslip(i)-1
        temp = height(seg_AE{i}{j});
        stack_FS_AE(k:k+temp-1,1:end-1) = seg_AE{i}{j};
        stack_FS_AE(k:k+temp-1,end) = i;
        stack_FS_time(k:k+temp-1) = kn_leverarm{i}(j);
        k = k+temp;
    end  
end

% FS duration measured based on the occurance of first FS AE
FS_duration = zeros(length(kn_leverarm),1);
for i = 1:length(farray)
    FS_duration(i) = kn_leverarm{i}(1);
end

% Summing the radii of all FS AE
FS_radius_sum = zeros(length(farray),1);
for i = 1:length(farray)
    temp = zeros(window_stickslip(i)-1,1);
    for j = 1:window_stickslip(i)-1
        %temp = height(seg_AE{i}{j});
        temp(j,1) = sum(seg_AE{i}{j}(:,7));
    end  
    FS_radius_sum(i,1) = FS_radius_sum(i,1)+sum(temp);
end
%
%stack_FS_loc, dip, strike, thickness, sqrt(dip^2+strike^2)
stack_FS_loc = (coord_trans * stack_FS_AE(:,4:6)')';
stack_FS = zeros(L,6);
stack_FS(:,1) = stack_FS_AE(:,end);%Simulation number
stack_FS(:,2) = stack_FS_AE(:,7);%AE radius
stack_FS(:,3) = stack_FS_AE(:,8);%AE magnitudes
stack_FS(:,4) = stack_FS_time;%time

% stacking mainshock
L = 0;
for i = 1:length(farray)
    L = L + height(seg_AE{i}{window_stickslip(i)});
end
W = length(seg_AE{1}{1}(1,:));
stack_MS_AE = zeros(L, W+1);
k = 1;
for i = 1:length(farray)
    temp = height(seg_AE{i}{window_stickslip(i)});
    stack_MS_AE(k:k+temp-1,1:end-1) = seg_AE{i}{window_stickslip(i)};
    stack_MS_AE(k:k+temp-1,end) = i;
    k = k+temp;
end
%
k = 1;
stack_MS_time = zeros(length(farray),1);
for i = 1:length(farray)
    temp = height(seg_AE{i}{window_stickslip(i)});
    %difference in time between each AE and the first AE
    AE_dt = abs(stack_MS_AE(k:k+temp-1,1)-stack_MS_AE(k,1));
    stack_MS_time(k:k+temp-1) = AE_dt;
    k = k+temp;
end

%time and location of the first MS AE
stack_MS_first = zeros(length(farray), 4);
k = 1;
for i = 1:length(farray)
    %MS AE time, x, y, z
    stack_MS_first(i,:) = stack_MS_AE(k,[2,4:6]);
    temp = height(seg_AE{i}{window_stickslip(i)});
    k = k+temp;
end

% Calculate the mainshock AE moment magnitude
seg_MS_moment = cell(length(farray),1);
seg_MS_moment_max = zeros(length(farray),1);
seg_MS_moment_sum = zeros(length(farray),1);
for i = 1:length(farray)
    %MS AE moment magnitude
    temp = seg_AE{i}{window_stickslip(i)};
    seg_MS_moment{i} = 10.^((temp(:,8)+6).*1.5);
    seg_MS_moment_max(i) = max(seg_MS_moment{i});
    seg_MS_moment_sum(i) = sum(seg_MS_moment{i});
end

%
stack_MS = zeros(L,6);
stack_MS(:,1) = stack_MS_AE(:,end);
stack_MS(:,2) = stack_MS_AE(:,7);
stack_MS(:,3) = stack_MS_AE(:,8); 
stack_MS(:,4) = stack_MS_time;

% For the distance between FS and the first MS AE, it requires 
% the knowledge of MS AE, and therefore, calculate the distances after MS 
% stacking

% How far each FS AE to the first MS AE
k = 1;
stack_FS_firstMS = zeros(height(stack_FS_AE),3);
for i = 1:length(farray)
    temp = sum(stack_FS_AE(:,19)==i);
    stack_FS_firstMS(k:k+temp-1,:) = stack_FS_loc(k:k+temp-1,1:3)...
        -(coord_trans*stack_MS_first(i,2:4)')';
    %distance
    stack_FS(k:k+temp-1,5) = ...
        vecnorm(stack_FS_firstMS(k:k+temp-1,:)')'...
        ;
    %angle (from 0 to 360 deg)
    stack_FS(k:k+temp-1,6) = rad2deg(...
        calculate_angle(...
            stack_FS_firstMS(k:k+temp-1,2),...
            stack_FS_firstMS(k:k+temp-1,1)...
            )...
        );
    k = k+temp;
end

% How far each MS AE to the first MS AE
k = 1;
for i = 1:length(farray)
    temp = height(stack_MS_AE(stack_MS_AE(:,19)==i,:));
    stack_MS_AE_each = stack_MS_AE(k:k+temp-1,:);
    %distance from subsequent AE to first
    stack_MS_firstMS = (coord_trans*stack_MS_AE_each(:,4:6)')'...
        -(coord_trans*stack_MS_first(i,2:4)')';
    stack_MS(k:k+temp-1,5) = vecnorm((stack_MS_firstMS'));
    stack_MS(k:k+temp-1,6) = rad2deg(...
        calculate_angle(...
            stack_MS_firstMS(:,1),stack_MS_firstMS(:,2)...
            )...
         );
    k = k+temp;
end

%
stack_FS_avg = zeros(length(farray), 5);
stack_FS_var = zeros(length(farray), 5);

for j = 2:6
    stack_FS_avg(1,j) = mean(stack_FS(stack_FS(:,1)==1, j));
    stack_FS_var(1,j) = var(stack_FS(stack_FS(:,1)==1, j));
    for i = 2:length(farray)
        stack_FS_avg(i,j) = (...
            mean(stack_FS(stack_FS(:,1)==i,j))+...
            stack_FS_avg(i-1,j)...
            )/2;
        stack_FS_var(i,j) = var( stack_FS(stack_FS(:,1)<=i, j) );
    end
end
%%
% plotting evolution of FS mean and var
fig = figure(f);f=f+1;
for i = 1:5
    subplot(2,3,i);
    plot(stack_FS_avg(:,i+1));
    set(gca,'FontName', 'Times');
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Increasing number of simulations');
ylabel('Evolution of means');
set(gca,'FontName', 'Times');

fig = figure(f);f=f+1;
for i = 1:5
    subplot(2,3,i);
    plot(stack_FS_var(:,i+1));
    set(gca,'FontName', 'Times');
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Increasing number of simulations');
ylabel('Evolution of variances');
set(gca,'FontName', 'Times');
%%

% Cross-correlation of FS variables
FS_corr = zeros(5,5,length(farray));
for k = 1:length(farray)
    for i = 2:6
        for j = i+1:6
            %FS_corrcov among every possible pairs of simulations
            A = stack_FS(stack_FS(:,1)==k,i);
            B = stack_FS(stack_FS(:,1)==k,j);
            temp = corrcov(cov(A,B));
            FS_corr(i-1,j-1,k) = temp(1,2);        
        end
    end
    %FS_corr was an upper triangle matrix
    %flip diagonal and add back to fill in the lower triangle
    FS_corr(:,:,k)=FS_corr(:,:,k)+flip(flip((FS_corr(:,:,k)'),1));
end
%
FS_corr_mean = zeros(5,5,length(farray));
FS_corr_var = zeros(5,5,length(farray));
for i = 1:length(farray)
    FS_corr_mean(:,:,i) = mean(FS_corr(:,:,1:i), 3);
    FS_corr_var(:,:,i) = var(FS_corr_mean,0,3);
end

%%
% Plot the evolution of the means of correlation coefficients of FS
% wrt to increasing number of simulations
fig=figure(f);f=f+1;
for i = 1:5
    for j = 1:5
        subplot(5,5,(i-1)*5+j);
        %flatten the 3rd dimension of FS_corr_mean for plot
        plt = reshape( FS_corr_mean(i,j,:),[],1 );
        plot(plt);
        set(gca,'FontName', 'Times');
    end
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Increasing number of simulations');
ylabel('Evolution of means of foreshock correlation coefficent');
set(gca,'FontName', 'Times');

% Plot the evolution of the variances of correlation coefficients of FS
% wrt to increasing number of simulations
fig=figure(f);f=f+1;
for i = 1:5
    for j = 1:5
        subplot(5,5,(i-1)*5+j);
        %flatten the 3rd dimension of FS_corr_mean for plot
        plt = reshape( FS_corr_var(i,j,:),[],1 );
        plot(plt);
        set(gca,'FontName', 'Times');
    end
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Increasing number of simulations');
ylabel('Evolution of variances of foreshock correlation coefficent');
set(gca,'FontName', 'Times');
%%
% Cross-correlation of MS variables
MS_corr = zeros(5,5,length(farray));
for k = 1:length(farray)
    for i = 2:6
        for j = i+1:6
            %MS_corrcov from each simulation
            A = stack_MS(stack_MS(:,1)==k,i);
            B = stack_MS(stack_MS(:,1)==k,j);
            temp = corrcov(cov(A,B));
            MS_corr(i-1,j-1,k) = temp(1,2);        
        end
    end
    %FS_corr was an upper triangle matrix
    %flip diagonal and add back to fill in the lower triangle
    MS_corr(:,:,k)=MS_corr(:,:,k)+flip(flip((MS_corr(:,:,k)'),1));
end

%
MS_corr_mean = zeros(5,5,length(farray));
MS_corr_var = zeros(5,5,length(farray));
for i = 1:length(farray)
    MS_corr_mean(:,:,i) = mean(MS_corr(:,:,1:i), 3);
    MS_corr_var(:,:,i) = var(MS_corr_mean,0,3);
end
%%
% Plot the evolution of the means of correlation coefficients of MS
% wrt to increasing number of simulations
fig=figure(f);f=f+1;
for i = 1:5
    for j = 1:5
        subplot(5,5,(i-1)*5+j);
        %flatten the 3rd dimension of FS_corr_mean for plot
        plt = reshape( MS_corr_mean(i,j,:),[],1 );
        plot(plt);
    end
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Increasing number of simulations');
ylabel('Evolution of means of mainshock correlation coefficent');
set(gca,'FontName', 'Times');

% Plot the evolution of the variances of correlation coefficients of MS
% wrt to increasing number of simulations
fig=figure(f);f=f+1;
for i = 1:5
    for j = 1:5
        subplot(5,5,(i-1)*5+j);
        %flatten the 3rd dimension of FS_corr_mean for plot
        plt = reshape( MS_corr_var(i,j,:),[],1 );
        plot(plt);
    end
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Increasing number of simulations');
ylabel('Evolution of variances of mainshock correlation coefficent');
set(gca,'FontName', 'Times');

%%

% Normalize the FS and MS vectors
stack_FS_norm = zeros(height(stack_FS), 6);
stack_FS_norm(:,1) = stack_FS(:,1);
stack_FS_norm(:,2) = normalize(stack_FS(:,2),'zscore');%AE radius
stack_FS_norm(:,3) = normalize(stack_FS(:,3),'zscore');%AE magnitudes
stack_FS_norm(:,4) = normalize(stack_FS(:,4),'zscore');%time
stack_FS_norm(:,5) = normalize(stack_FS(:,5),'zscore');%dist to first MS AE
stack_FS_norm(:,6) = normalize(stack_FS(:,6),'zscore');%direction to first MS AE

stack_MS_norm = zeros(height(stack_MS),6);
stack_MS_norm(:,1) = stack_MS_AE(:,end);
stack_MS_norm(:,2) = normalize(stack_MS(:,2),'zscore');
stack_MS_norm(:,3) = normalize(stack_MS(:,3),'zscore'); 
stack_MS_norm(:,4) = normalize(stack_MS(:,4),'zscore');
stack_MS_norm(:,5) = normalize(stack_MS(:,5),'zscore'); 
stack_MS_norm(:,6) = normalize(stack_MS(:,6),'zscore');

%% Calculate Hausdorff distances among FS and MS, 
% and apply hierarchy clustering
n = 2;
H = zeros(length(farray), length(farray));%upper diagnal of distance matrix
for j = 1:length(farray)
    for i = j:length(farray)
        H(j,i) = ...
            hausdorff( stack_FS(stack_FS(:,1)==j,n:end),...
                stack_FS(stack_FS(:,1)==i,n:end),...
                inv(cov(stack_FS(:,n:end))) ); 
    end
end
FS_hausdorff = H + H' - diag([diag(H)]);
fsZ = linkage(FS_hausdorff,'complete');
Z = figure(f);f=f+1;
dendrogram(fsZ);
xlabel('Simulation ID');
ylabel('Hausdorff distance');
title('Clustering based on foreshocks using complete linkage');
set(gca,'FontName', 'Times');
yticks( round(min(fsZ(:,end))/2)*2:2 ...
    :round(max(fsZ(:,end))/2)*2 );
figurename = strcat(fgroup,'FSlinkage_HD_COM.jpg'); 
%exportgraphics(Z,figurename,'Resolution',300);%saveas(figure1,figurename);

H = zeros(length(farray), length(farray));%upper diagnal of distance matrix
for j = 1:length(farray)
    for i = j:length(farray)
        H(j,i) = ...
            hausdorff( stack_MS(stack_MS(:,1)==j,n:6),...
            stack_MS(stack_MS(:,1)==i,n:6), ...
            inv(cov(stack_MS(:,n:end))) ); 
    end
end
MS_hausdorff = H + H' - diag([diag(H)]);
msZ = linkage(MS_hausdorff,'complete'); 
Z = figure(f);f=f+1;
dendrogram(msZ);
xlabel('Simulation ID');
ylabel('Hausdorff distance');
title('Clustering based on mainshocks using complete linkage');
set(gca,'FontName', 'Times');
yticks( round(min(msZ(:,end))/2)*2:2 ...
    :round(max(msZ(:,end))/2)*2 );
figurename = strcat(fgroup,'MSlinkage_HD_COM.jpg'); %saveas(figure1,figurename);
%exportgraphics(Z,figurename,'Resolution',300);

% Check, with increasing number of simulations, whether average distance
% among any pair of simulations (resemblance) changes
FS_hausdorff_evol_mean = zeros(length(farray),length(farray)-1);
for i = 2:length(farray)
    study_data = FS_hausdorff(:, 1:i);
    FS_hausdorff_evol_mean(:,i-1)=mean(study_data, 2)/(i-1)*i;
end
%
fig = figure(f); f=f+1;
for i = 1:length(farray)
    subplot(5,6,i);
    plot(2:length(farray),FS_hausdorff_evol_mean(i,:));

end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Increasing number of simulations');
ylabel({'Evolution of the mean pair-wise Hausdorff distance'});
title('Foreshock');
set(gca,'FontName', 'Times');
figurename = strcat(fgroup,'FS_evol_mean_HD.jpg'); 
exportgraphics(fig,figurename,'Resolution',300);
%
% Calculate FM Index among the two hierarchical trees
cluster_eval = zeros(length(farray)-1,1);
m_k = cell(length(farray),1);
for i = 2:length(farray)
    a = cluster(fsZ, 'maxclust', i);
    b = cluster(msZ, 'maxclust', i);
    [cluster_eval(i-1),m_k{i-1}] = Fowlkes_Mallows_Index(a, b);
end
% figure(f);f=f+1;
% plot(cluster_eval);
% xlabel('Increasing number of clusters');
% ylabel('Fowlkes Mallows Index');
% set(gca,'FontName', 'Times');

% %%
% %scatter(log10(stack_FS_AE(:,7)), stack_FS_AE(:,8));
% for i = 1:2%length(farray)
%     x = stack_FS_AE(stack_FS_AE(:,end)==i,7);
%     y = stack_FS_AE(stack_FS_AE(:,end)==i,8);
%     z = stack_FS_time(stack_FS_AE(:,end)==i);
%     scatter3(x, y, z);
%     hold on
% end
% hold off

% FM Index B(k,k) p test
Bk_null = zeros(length(farray),3);
for i = 2:length(farray)-1 %from k = 2 to k = 26
    [Bk_null(i,1),Bk_null(i,2),Bk_null(i,3)] ...
        = FM_null_hypothesis_test(i, fsZ, msZ);
end


%% Analysis based on Chamfer distance
CD = zeros(length(farray), length(farray));%upper diagnal of distance matrix
for j = 1:length(farray)
    for i = j:length(farray)
        a = stack_FS_norm(stack_FS_norm(:,1)==j,n:end);
        b = stack_FS_norm(stack_FS_norm(:,1)==i,n:end);
        
        pdist_ab_ab = squareform(pdist([a;b]));
        pdist_a_b = pdist_ab_ab(1:height(a), height(a)+1:end);
        
        CD(j,i) = (min(pdist_a_b,[],2))'*(min(pdist_a_b,[],2)) + ...
            (min(pdist_a_b,[],1))*(min(pdist_a_b,[],1))';
    end
end
FS_chamfer = CD + CD' - diag([diag(CD)]);
fsZ_CD = linkage(FS_chamfer, 'complete');
Z = figure(f);f=f+1;
dendrogram(fsZ_CD);
xlabel('Simulation ID');
ylabel('Chamfer distance');
title('Clustering based on foreshocks using complete linkage');
set(gca,'FontName', 'Times');
yticks(round(min(fsZ_CD(:,end))/100)*100:...
    100:...
    round(max(fsZ_CD(:,end))/100)*100 );
figurename = strcat(fgroup,'FSlinkage_CD_COM.jpg'); 
%exportgraphics(Z,figurename,'Resolution',300);
%
CD = zeros(length(farray), length(farray));%upper diagnal of distance matrix
for j = 1:length(farray)
    for i = j:length(farray)
        a = stack_MS_norm(stack_MS_norm(:,1)==j,n:end);
        b = stack_MS_norm(stack_MS_norm(:,1)==i,n:end);
        
        pdist_ab_ab = squareform(pdist([a;b]));
        pdist_a_b = pdist_ab_ab(1:height(a), height(a)+1:end);
        
        CD(j,i) = (min(pdist_a_b,[],2))'*(min(pdist_a_b,[],2)) + ...
            (min(pdist_a_b,[],1))*(min(pdist_a_b,[],1))';
    end
end
MS_chamfer = CD + CD' - diag([diag(CD)]);
msZ_CD = linkage(MS_chamfer, 'complete');
Z = figure(f);f=f+1;
dendrogram(msZ_CD);
xlabel('Simulation ID');
ylabel('Chamfer distance');
title('Clustering based on mainshocks using complete linkage');
set(gca,'FontName', 'Times');
yticks(round(min(msZ_CD(:,end))/100)*100:...
    100:...
    round(max(msZ_CD(:,end))/100)*100 );
figurename = strcat(fgroup,'MSlinkage_CD_COM.jpg'); 
%exportgraphics(Z,figurename,'Resolution',300);

% Check, with increasing number of simulations, whether average distance
% among any pair of simulations (resemblance) changes
FS_Chamfer_evol_mean = zeros(length(farray),length(farray)-1);
for i = 2:length(farray)
    study_data = FS_chamfer(:, 1:i);
    FS_Chamfer_evol_mean(:,i-1)=mean(study_data, 2)/(i-1)*i;
end
%
fig = figure(f); f=f+1;
for i = 1:length(farray)
    subplot(5,6,i);
    plot(2:length(farray),FS_Chamfer_evol_mean(i,:));

end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Increasing number of simulations');
ylabel({'Evolution of the mean pair-wise Chamfer distance'});
title('Foreshock');
set(gca,'FontName', 'Times');
figurename = strcat(fgroup,'FS_evol_mean_CD.jpg'); 
exportgraphics(fig,figurename,'Resolution',300);

% Calculate FM Index among the two hierarchical trees
cluster_eval_CD = zeros(length(farray)-1,1);
m_k_CD = cell(length(farray),1);
for i = 2:length(farray)
    a = cluster(fsZ_CD, 'maxclust', i);
    b = cluster(msZ_CD, 'maxclust', i);
    [cluster_eval_CD(i-1),m_k_CD{i-1}] = Fowlkes_Mallows_Index(a, b);
end
% FM Index B(k,k) p test
Bk_null_CD = zeros(length(farray),3);
for i = 2:length(farray)-1 %from k = 2 to k = 26
    [Bk_null_CD(i,1),Bk_null_CD(i,2),Bk_null_CD(i,3)] ...
        = FM_null_hypothesis_test(i, fsZ_CD, msZ_CD);
end

%% Earth Mover's distance
EMD = zeros(length(farray), length(farray));
for j = 1:length(farray)
    for i = j:length(farray)
        a = stack_FS_norm(stack_FS_norm(:,1)==j,n:end);
        b = stack_FS_norm(stack_FS_norm(:,1)==i,n:end);
        
        if height(a) < height(b)
            c = b;
            b = a;
            a = c;
        end
        
        pdist_ab_ab = squareform(pdist([a;b]));
        pdist_a_b = pdist_ab_ab(1:height(a), height(a)+1:end);

        min_dist = zeros(size(pdist_a_b,2),1);
        for ii = 1:size(pdist_a_b,2)
            min_dist(ii) = min(min(pdist_a_b));
            [min_row, min_col] = find(pdist_a_b == min_dist(ii),1,'first');
            pdist_a_b(min_row,:) = [];
            pdist_a_b(:,min_col) = [];
        end
        EMD(j,i) = sum(min_dist)/length(b);
    end
end
FS_EMD = EMD + EMD' - diag([diag(EMD)]);
fsZ_EMD = linkage(FS_EMD, 'complete');
Z = figure(f);f=f+1;
dendrogram(fsZ_EMD);
xlabel('Simulation ID');
ylabel('Earth mover distance');
title('Clustering based on foreshocks using complete linkage');
set(gca,'FontName', 'Times');
yticks( round(min(fsZ_EMD(:,end))*2)/2 : 0.5 : ...
    round(max(fsZ_EMD(:,end))*2)/2 );
figurename = strcat(fgroup,'FSlinkage_EMD_COM_V2.jpg'); 
% exportgraphics(Z,figurename,'Resolution',300);

EMD = zeros(length(farray), length(farray));
for j = 1:length(farray)
    for i = j:length(farray)
        a = stack_MS_norm(stack_MS_norm(:,1)==j,n:end);
        b = stack_MS_norm(stack_MS_norm(:,1)==i,n:end);
        
        if height(a) < height(b)
            c = b;
            b = a;
            a = c;
        end
        
        pdist_ab_ab = squareform(pdist([a;b]));
        pdist_a_b = pdist_ab_ab(1:height(a), height(a)+1:end);

        min_dist = zeros(size(pdist_a_b,2),1);
        for ii = 1:size(pdist_a_b,2)
            min_dist(ii) = min(min(pdist_a_b));
            [min_row, min_col] = find(pdist_a_b == min_dist(ii),1,'first');
            pdist_a_b(min_row,:) = [];
            pdist_a_b(:,min_col) = [];
        end
        EMD(j,i) = sum(min_dist)/length(b);
    end
end
MS_EMD = EMD + EMD' - diag([diag(EMD)]);
msZ_EMD = linkage(MS_EMD, 'complete');
Z = figure(f);f=f+1;
dendrogram(msZ_EMD);
xlabel('Simulation ID');
ylabel('Earth mover distance');
title('Clustering based on mainshocks using complete linkage');
set(gca,'FontName', 'Times');
yticks( round(min(msZ_EMD(:,end))*2)/2 : 0.5 : ...
    round(max(msZ_EMD(:,end))*2)/2 );
figurename = strcat(fgroup,'MSlinkage_EMD_COM_V2.jpg'); 
% exportgraphics(Z,figurename,'Resolution',300);

%% 
% Check, with increasing number of simulations, whether average distance
% among any pair of simulations (resemblance) changes
FS_EMD_evol_mean = zeros(length(farray),length(farray)-1);
for i = 2:length(farray)
    study_data = FS_EMD(:, 1:i);
    FS_EMD_evol_mean(:,i-1)=mean(study_data, 2)/(i-1)*i;
end
%
fig = figure(f); f=f+1;
for i = 1:length(farray)
    subplot(5,6,i);
    plot(2:length(farray),FS_EMD_evol_mean(i,:));

end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('Increasing number of simulations');
ylabel({'Evolution of the mean pair-wise Earth mover distance'});
title('Foreshock');
set(gca,'FontName', 'Times');
figurename = strcat(fgroup,'FS_evol_mean_EMD.jpg'); 
exportgraphics(fig,figurename,'Resolution',300);
%%
% Calculate FM Index among the two hierarchical trees
cluster_eval_EMD = zeros(length(farray)-1,1);
m_k_EMD = cell(length(farray),1);
for i = 2:length(farray)
    a = cluster(fsZ_EMD, 'maxclust', i);
    b = cluster(msZ_EMD, 'maxclust', i);
    [cluster_eval_EMD(i-1),m_k_EMD{i-1}] = Fowlkes_Mallows_Index(a, b);
end
% FM Index B(k,k) p test
Bk_null_EMD = zeros(length(farray),3);
for i = 2:length(farray)-1 %from k = 2 to k = 26
    [Bk_null_EMD(i,1),Bk_null_EMD(i,2),Bk_null_EMD(i,3)] ...
        = FM_null_hypothesis_test(i, fsZ_EMD, msZ_EMD);
end
%%
%Plot distribution of data
fig = figure(f); f=f+1;
plt_label = {'AE radius [mm]','AE magnitudes',...
    'Time [steps]','Distance [mm]','Direction [deg]'};
num_bins = 20; % Specify the number of bins

% Define the figure dimensions with a ratio of 1:5
figure_width = 2000;
figure_height = 400;
set(fig, 'Position', [100, 100, figure_width, figure_height]);

for i = 1:5
    subplot(1,5,i);
    if i == 3 
        histogram(-stack_FS(:,i+1), 'NumBins', num_bins);
    else
        histogram(stack_FS(:,i+1), 'NumBins', num_bins);
    end
    xlabel(plt_label(i));
    set(gca,'FontSize',12,'FontName','Times');
    ylim([0,700]);
end

exportgraphics(fig,'FS_distribution_231022.jpg','Resolution',300);

fig=figure(f); f=f+1;
plt_label = {'AE radius [mm]','AE magnitudes',...
    'Time [steps]','Distance [mm]','Direction [deg]'};
set(fig, 'Position', [100, 100, figure_width, figure_height]);

for i = 1:5
    subplot(1,5,i);
    histogram(stack_MS(:,i+1), 'NumBins', num_bins);
    xlabel(plt_label(i));
    set(gca,'FontSize', 12,'FontName', 'Times');
    ylim([0,450]);
end

exportgraphics(fig,'MS_distribution_231022.jpg','Resolution',300);
%
% PCA
% [fscoeff,fsscore,fslatent,fstsquared,fsexplained]...
%     = pca(stack_FS_norm(:,2:6));
%% Grain-scale test
%checking to see if the jump ladders and other grainscale observations
%show any temporal evolution trends in their means and variances
S6_grainscale
%% Evolution of means
S6_evolutionofmeans

%% cluster based on shape of the CN segments
S6_stackingCNshape

%% Temporal evolution of particle speed
S6_tempEvolParticleSpeed

%% cluster each stack_JL_FS variables using kmeans
Stackingdata_V6_cluster_kmean

%% Plot foreshock and mainshock AE from Simulation 6 and 7
S6_scatterFSMSAE