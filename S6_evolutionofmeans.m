%Stackingdata_v6 subscript
% Evolution of means
stack_JL_FS_avg = zeros(length(stack_JL_FS),9);
stack_JL_FS_var = zeros(length(stack_JL_FS),9);
for j = 1:9
    if j == 4 || j==5, continue;  end
    avg = 0;
    new = 0;
    sum_sq = 0;
    for i = 1:length(stack_JL_FS)
        new = stack_JL_FS(i,j);
        avg = [(i-1)/i, 1/i]*[avg; new];
        stack_JL_FS_avg(i,j) = avg;
        
        sum_sq = new^2 + sum_sq;
        new_var = sum_sq/i - avg^2;
        stack_JL_FS_var(i,j) = new_var;
    end
    figure(f);f=f+1;plot(stack_JL_FS_var(:,j));
end
%
stack_JL_bcoord_avg = zeros(length(stack_JL_bcoord),3);
stack_JL_bcoord_var = zeros(length(stack_JL_bcoord),3);
for j = 1:3
    
    avg = 0;
    new = 0;
    sum_sq = 0;
    for i = 1:length(stack_JL_FS)
        new = stack_JL_FS(i,j+3);
        avg = [(i-1)/i, 1/i]*[avg; new];
        stack_JL_bcoord_avg(i,j) = avg;
        
        sum_sq = new^2 + sum_sq;
        new_var = sum_sq/i - avg^2;
        stack_JL_bcoord_var(i,j) = new_var;
    end
    figure(f);f=f+1;plot(stack_JL_bcoord_var(:,j));
end
%
for j = 1:9
    if j == 4 || j==5, continue;  end
    
    figure(f);f=f+1;plot(stack_JL_FS_avg(:,j));
end