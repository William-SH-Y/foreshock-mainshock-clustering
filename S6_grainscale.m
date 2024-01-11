%stackigndata_V6 subscript
%grain-scale test
%checking to see if the jump ladders and other grainscale observations
%show any temporal evolution trends in their means and variances

k = 1;
L = 0;
for j = 1:length(farray) %looping through simulations
    for i = 1:window_stickslip(j)-1 %looping through foreshock episodes
        stack_height = size(Bcoord{j},1);
        L = L + stack_height;
    end
end
%
stack_JL_bcoord = zeros(L,9);
k = 1;
% concatenate simulation ID, episode ID, sample ID and bcoord
for j = 1:length(farray) %looping through simulations
    for i = 1:window_stickslip(j)-1 %looping through foreshock episodes
        stack_height = size(Bcoord{j},1);
        stack_JL_bcoord(k:k+stack_height-1,1:3) = ...
            [ones(stack_height,1)*j,... %sim ID
             ones(stack_height,1)*i,... %episode ID
             (1:1:stack_height)',... %sample ID
            ]; 
        stack_JL_bcoord(k:k+stack_height-1,4:end) = ...
            Bcoord{j};
        k = k+stack_height;
    end
end

k = 1;
stack_JL_FS = zeros(L,15);
for j = 1:length(farray) %looping through simulations
    for i = 1:window_stickslip(j)-1 %looping through foreshock episodes
        stack_height = size(Bcoord{j},1);
        stack_JL_FS(k:k+stack_height-1,:) = ...
            [ones(stack_height,1)*kn_leverarm{j}(i), ... %time
             jump_Ladder(j).cn(:,i),...
             jump_Ladder(j).cs(:,i),...
             jump_Cat(j).cn(:,i),...
             jump_Cat(j).cs(:,i),...
             jump_Rat(j).prejump_cn(:,i),...
             jump_Rat(j).prejump_cs(:,i),...
             jump_Rat(j).postjump_cn(:,i),...
             jump_Rat(j).postjump_cs(:,i),...
             transpose(jump_Coef(j).prejump_bdx(i,:,1)),... %10,prejump bdx
             transpose(jump_Coef(j).prejump_bdy(i,:,1)),...
             transpose(jump_Coef(j).prejump_bdz(i,:,1)),...
             transpose(jump_Coef(j).postjump_bdx(i,:,1)),...
             transpose(jump_Coef(j).postjump_bdy(i,:,1)),...
             transpose(jump_Coef(j).postjump_bdz(i,:,1))];
        k = k+stack_height;
    end
end

% remove nan rows
idx_nan = find(isnan(stack_JL_FS) == 1);
[idx_nan_row,idx_nan_col] = ind2sub(size(stack_JL_FS),idx_nan);
if sum(sum(isnan(stack_JL_FS)==1))>0 %only execute if there is nan
    stack_JL_FS(idx_nan_row,:)=[];
    stack_JL_bcoord(idx_nan_row,:)=[];
end
% remove inf rows
idx_inf = find(isinf(stack_JL_FS) == 1);
[idx_inf_row,idx_inf_col] = ind2sub(size(stack_JL_FS),idx_inf);
if sum(sum(isinf(stack_JL_FS)==1))>0 %only execute if there is nan
    stack_JL_FS(idx_inf_row,:)=[];
    stack_JL_bcoord(idx_inf_row,:)=[];
end

% remove outliers
[stack_JL_FS, TF] = rmoutliers(stack_JL_FS, 'percentiles', [1,99]);
stack_JL_bcoord = stack_JL_bcoord(~TF,:);

% normalize
stack_JL_FS_norm = zeros(size(stack_JL_FS));
for i = 1:size(stack_JL_FS,2)
    if i == 4 || i == 5
        stack_JL_FS_norm(:,i) = stack_JL_FS(:,i); 
    else
        stack_JL_FS_norm(:,i) = normalize(stack_JL_FS(:,i),'zscore');
    end
end

% setting up time segments for temporal evolution calculation
t = linspace(max(stack_JL_FS(:,1)),min(stack_JL_FS(:,1)),10);
stack_JL_FS_timeseg = zeros(height(stack_JL_FS), 1);
for j = 2:length(t)
    cond = stack_JL_FS(:,1)>=t(j) & ...
            stack_JL_FS(:,1)<=t(j-1); %within the time segment
    stack_JL_FS_timeseg(cond) = j-1;
end