

clear
close all
clc
% importing
f = 1;
mv_W = 40e-3;
mv_H = 107e-3;

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
    
    load(strcat('S',farray{i},'broken_trigger'),'broken_trigger');
    broken_Trigger{i} = broken_trigger; %stickslip stressdrop
    
    load(strcat('S',farray{i},'broken_stresstrigger'),'broken_stresstrigger');
    broken_ST{i} = broken_stresstrigger; %stickslip stressdrop
    load(strcat('S',farray{i},'broken_inttriggerx'),'broken_inttriggerx');
    broken_IXT{i} = broken_inttriggerx; %stickslip stressdrop
    load(strcat('S',farray{i},'broken_inttriggery'),'broken_inttriggery');
    broken_IYT{i} = broken_inttriggery; %stickslip stressdrop
    load(strcat('S',farray{i},'broken_inttriggerz'),'broken_inttriggerz');
    broken_IZT{i} = broken_inttriggerz; %stickslip stressdrop
    
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

for i = 1:length(farray)
    temp = zeros(size(seg_Sjcrk{i},1), 1);
    for j = 1:size(seg_Sjcrk{i},1)
        temp(j) = size(seg_Sjcrk{i}{j},1);
    end
    window_stickslip(i) = find(temp == max(temp));
    clear temp
end
%% Stacking all the particle velocities (from raw bd data) such that the trend of acceleration can be seen
bd_stack = cell(length(farray),1);
bd_stackweighted = cell(length(farray),1);
bd_weight = zeros(length(farray),6);
for i = 1:length(farray) 
    temp = (1:find(cyc_Wrtss{i} == 0));
    temp1 = bd_Mean{i}(temp, :);
    bd_stack{i} = flip(temp1,1);
    bd_weight(i,:) = bd_Size(i,:) ./ sum(bd_Size);
    
    for j = 1:6,bd_stackweighted{i}(:,j)=bd_stack{i}(:,j)*bd_weight(i,j);end
    

end
%
bd_stacklength = ...
    min(cell2mat(cellfun(@length, bd_stackweighted, 'UniformOutput', false)));
bd_stackweightedtruncate = zeros(bd_stacklength,6,length(farray));
for i = 1:length(farray)
    bd_stackweightedtruncate(:,:,i) = bd_stackweighted{i}(1:bd_stacklength,:);
end
bd_stackresult = sum(bd_stackweightedtruncate, 3);
%
figure(f);f=f+1;
for i = 1:6
    subplot(2,3,i);
    plot(flip(bd_stackresult(:,i)));
end

%% FS
% cross covariance matrix
%remove the NaN rows
shock = cell(length(farray),1);
FS = cell(length(farray),1);
for i = 1:length(farray)
    shock{i} = [kn_leverarm{i},kn_Pks{i}',sj_spacdist{i}];
    shock{i}(isnan(shock{i}(:,3)), : ) = [];
    FS{i} = shock{i}((1:find(shock{i}(:,1)==0)-1 ),: );
end

for i = 1:length(farray)
    FSbar(i,:) = mean(FS{i},1);
    FSvar(i,:) = var(FS{i},1);
    FScov{i} = ((FS{i} - FSbar(i,:))' * (FS{i} - FSbar(i,:)))/ (length(FS{i}(:,1))-1); %same as cov(FS{i})    
    [FScorr{i},s(i,:)] = corrcov(FScov{i}); %Use matlab built in function to calculate correlation matrix
end
%
%writematrix(FScorr{end},strcat(fgroup,'FScorr.csv'));

% calculate average and variance of FScorr with respect to number of
% simulations (to see if it converges)
for ii = 1:length(farray)
    for i = 1:size(FS{1},2) %1:5, 5 dimensions
        for j = 1:size(FS{1},2) %1:5, 5 dimensions
            for k = 1:ii%length(farray)
                temp(k) = FScorr{k}(i,j);
            end
            FScorravg(i,j, ii) = mean(temp);
            FScorrvar(i,j, ii) = var(temp);
        end
    end
end

figure1 = figure(f); f=f+1;
for i = 1:size(FS{1},2)
    for j = 1:size(FS{1},2)
        temp = zeros(length(farray),1);
        subplot(size(FS{1},2), size(FS{1},2), (size(FS{1},2)*i-(size(FS{1},2)-j)) );
        temp(:,1) = FScorravg(i, j, :);
        plot(temp);
        clear temp
    end
end


% Find the prominant peaks in foreshock kernel density function. First
% reorder the peaks, and then pick the top 5
FSreorder = cell(length(farray),1);
for i = 1:length(farray)
    FSreorder{i} = sortrows(FS{i}, 2, 'descend');
end
% Compile the prominant peaks in FS kernel density function
tempn = 5;
FShigh = zeros(tempn, 5, length(farray));
for i = 1:length(farray)
    FShigh(:,:,i) = FSreorder{i}(1:tempn,:);
end
FShightotal = zeros(tempn*length(farray),5);
j = 1;
for i = 1:length(farray)
    FShightotal(j:(j+tempn-1), 1:5) = FShigh(:,:,i);
    FShightotal(j:(j+tempn-1), 6) = i;
    j = j+tempn;
end

% Calculate the total mean and variance of all foreshocks
Ltotal = size(FS{1},1);
for i = 2:length(farray)
    Ltotal = Ltotal + size(FS{i},1);
end
FStotal = zeros( Ltotal, 6);
pltl = size(FStotal,1);
j = 1;
for i = 1:length(farray)
    FStotal(j:j+length(FS{i}(:,1))-1, 1:5) = FS{i};
    FStotal(j:j+length(FS{i}(:,1))-1, 6) = i;
    FStotal(j:j+length(FS{i}(:,1))-1, 7) = 1:size(FS{i},1);
    j = j+length(FS{i}(:,1));
end
FStotalmean = mean(FStotal,1);
FStotalvar = var(FStotal,1);

FStotalcov = cov(normalize(FStotal(:,1:5), 'range')');
FStotalcor = corrcov(FStotalcov);




%% (keyword stacking) Extract similar foreshocks from FShigh & Data stacking
j = 1;
for i = 1:length(farray)
    FShightotal(j:(j+tempn-1), 7) = 1:tempn;
    j = j+tempn;
end
FShtidvcor = corrcov(corr ( normalize(FShightotal(:,1:5),'range')' ) );
%
pltl = size(FShightotal,1);
pltx = linspace(-ceil(pltl/2), floor(pltl/2), pltl);
[pltX, pltY] = meshgrid(pltx, pltx);
FShtidvcorh = FShtidvcor;
FShtidvcorh(FShtidvcorh < 0.995) = 0;
%FShtidvcorh(FShtidvcorh < 0.99) = 0;
%FShtidvcorh(FShtidvcorh > 0.01 | FShtidvcorh < -0.01) = 0;
%figure(f); f=f+1; surf(pltX, pltY, FShtidvcorh); 
%
clearvars temprow tempcol temp1 temp2 tempi FS_snpn1 FS_snpn2
[temprow, tempcol] = find(FShtidvcorh);
temp1 = zeros(size(temprow,1),2);
temp2 = zeros(size(tempcol,1),2);
temp1(:,1) = FShightotal(temprow,6); temp1(:,2) = FShightotal(temprow,7);
temp2(:,1) = FShightotal(tempcol,6); temp2(:,2) = FShightotal(tempcol,7);

% Remove the diagnol from correlation matrix
temp = zeros(size(temp1,1),1);
for i = 1:size(temp1,1)
    temp(i) = isequal([temp1(i,:)], [temp2(i,:)]);
end
tempi = find(temp);
temp1(tempi,:) = []; temp2(tempi,:) = []; % Remove the diagnol from correlation matrix
%
clear tempi
tempi = zeros(size(temp1,1),2); % Redefine tempi as a temporary index 
for i = 1:size(temp1,1)
    tempi(i,1) = find(FShightotal(:,6) == temp1(i,1) & FShightotal(:,7) == temp1(i,2));
    tempi(i,2) = find(FShightotal(:,6) == temp2(i,1) & FShightotal(:,7) == temp2(i,2));
end
%
FS_simn_pkn = zeros(size(temp1,1), 5);
for i = 1:size(temp1,1)
    tempgroup = FShightotal(tempi(i,1), 6);
    tempc = FShightotal(tempi(i,1), 1);
    temptemp = find(FS{tempgroup}(:,1) == tempc);
    FS_simn_pkn(i,1:2) = [tempgroup, temptemp];
    
    tempgroup = FShightotal(tempi(i,2), 6);
    tempc = FShightotal(tempi(i,2), 1);
    temptemp = find(FS{tempgroup}(:,1) == tempc);
    FS_simn_pkn(i,3:4) = [tempgroup, temptemp];
end
% Separate the foreshocks based on x location
for i = 1:size(temp1, 1)
    tempx = FShightotal(tempi(i,1), 3);
    if tempx >= 0
        FS_simn_pkn(i,5) = 1;
    else
        FS_simn_pkn(i,5) = 0;
    end
end
%
FS_snpn1 = FS_simn_pkn(FS_simn_pkn(:,5)==1,:); %foreshock simulation number KDF peak number
FS_snpn2 = FS_simn_pkn(FS_simn_pkn(:,5)==0,:);

%
FStestgroupidx =  find(( ...
     FStotal(:,3)<-0.002 ...
)); 
FS_snpn1 = FStotal(FStestgroupidx, 6:7);
FStestgroupidx =  find(( ...
     FStotal(:,3)>0.002 ...
)); 
FS_snpn2 = FStotal(FStestgroupidx, 6:7);

%FStestgroup = FStotal(FStotal(:,1)>2e5 & FStotal(:,1)<3e5 & FStotal(:,2)<0.5e-6 & FStotal(:,3)>0, :);
FStestgroupidx =  find(( FStotal(:,1)<2e5 &...
     FStotal(:,3)<0.001 & FStotal(:,3)>-0.001 ...
)); 
% FStestgroupidx = find(FStotal(:,1)<5e5);    & FStotal(:,4)<0.005 &
% FStotal(:,4)>-0.005 FStotal(:,2)<1e-6 & 
FS_snpn3 = FStotal(FStestgroupidx, 6:7);

FS_snpn0 = [[1:length(farray)]' , window_stickslip'];



%% keyword stacking: part 1.5
ind1_0=zeros(size(FS_snpn0,1),1); ind2_0=zeros(size(FS_snpn0,1),1); leng_0=zeros(size(FS_snpn0,1),1);
for i = 1:size(FS_snpn0)
    sn = FS_snpn0(i,1);
    leng_0(i) = size(Bcoord{sn},1);
    if i == 1 
        ind1_0(i) = 1+0; ind2_0(i) = ind1_0(i)+leng_0(i)-1;
    else
        ind1_0(i) = 1+ind2_0(i-1); ind2_0(i) = ind1_0(i)+leng_0(i)-1;
    end     
end
mem_alloc_rows_0 = ind2_0(end);
%
stack_jl0 = zeros(mem_alloc_rows_0, length(fn)); %stacking jump ladder, fn:field names, fields ranging from bdx, dby, etc... so a total of 11
stack_jc0 = zeros(mem_alloc_rows_0, length(fn)); %stacking jump category
stack_jr0 = zeros(mem_alloc_rows_0, length(fnR)); %stacking jump ratio, fn:field names for jump Ratio
stack_bc0 = zeros(mem_alloc_rows_0, 6); %stacking bcoord
stack_cc0 = zeros(mem_alloc_rows_0, 3); %stacking ccoord
stack_Npre0 = zeros(mem_alloc_rows_0, 2);
% stack_broke1 = zeros(mem_alloc_rows_1, 20); %stacking ccoord
for i = 1:size(FS_snpn0,1)
    clear temp
    temp = Bcoord{FS_snpn0(i,1)};
    stack_bc0(ind1_0(i):ind2_0(i), :) = temp;
    clear temp
    for j = 1:length(fn)
        temp = jump_Ladder(FS_snpn0(i,1)).(fn{j})(:,FS_snpn0(i,2)) ;
        stack_jl0(ind1_0(i):ind2_0(i), j) = temp;
    end 
    clear temp
    for j = 1:length(fnR)
        temp = jump_Rat(FS_snpn0(i,1)).(fnR{j})(:,FS_snpn0(i,2)) ;
        stack_jr0(ind1_0(i):ind2_0(i), j) = temp;
    end 
    clear temp
    for j = 1:length(fn)
        temp = jump_Cat(FS_snpn0(i,1)).(fn{j})(:,FS_snpn0(i,2)) ;
        stack_jc0(ind1_0(i):ind2_0(i), j) = temp;
    end  
    
    stack_Npre0(ind1_0(i):ind2_0(i), 1) = seg_CN_pre{FS_snpn0(i,1)}(:,FS_snpn0(i,2));
    stack_Npre0(ind1_0(i):ind2_0(i), 2) = seg_CS_pre{FS_snpn0(i,1)}(:,FS_snpn0(i,2));
%     
%     clear temp
%     for j = 1:20
%         temp = broken_Trigger{FS_snpn1(i,1)}{FS_snpn1(i,2),j} ;
%         stack_broke1(ind1_1(i):ind2_1(i), j) = temp;
%     end 
    
end
stack_cc0 = [(stack_bc0(:,1)+stack_bc0(:,4))/2 , ...
    (stack_bc0(:,2)+stack_bc0(:,5))/2, (stack_bc0(:,3)+stack_bc0(:,6))/2];
stack_jr0(isnan(stack_jr0)) = 0;
stack_jr0(isinf(stack_jr0)) = 0;
%

ind1_1=zeros(size(FS_snpn1,1),1); ind2_1=zeros(size(FS_snpn1,1),1); leng_1=zeros(size(FS_snpn1,1),1);
for i = 1:size(FS_snpn1)
    sn = FS_snpn1(i,1);
    leng_1(i) = size(Bcoord{sn},1);
    if i == 1 
        ind1_1(i) = 1+0; ind2_1(i) = ind1_1(i)+leng_1(i)-1;
    else
        ind1_1(i) = 1+ind2_1(i-1); ind2_1(i) = ind1_1(i)+leng_1(i)-1;
    end     
end
mem_alloc_rows_1 = ind2_1(end);
%
stack_jl1 = zeros(mem_alloc_rows_1, length(fn)); %stacking jump ladder, fn:field names, fields ranging from bdx, dby, etc... so a total of 11
stack_jc1 = zeros(mem_alloc_rows_1, length(fn)); %stacking jump category
stack_jr1 = zeros(mem_alloc_rows_1, length(fnR)); %stacking jump ratio, fn:field names for jump Ratio
stack_bc1 = zeros(mem_alloc_rows_1, 6); %stacking bcoord
stack_cc1 = zeros(mem_alloc_rows_1, 3); %stacking ccoord
stack_Npre1 = zeros(mem_alloc_rows_1, 2);
% stack_broke1 = zeros(mem_alloc_rows_1, 20); %stacking ccoord
for i = 1:size(FS_snpn1,1)
    clear temp
    temp = Bcoord{FS_snpn1(i,1)};
    stack_bc1(ind1_1(i):ind2_1(i), :) = temp;
    clear temp
    for j = 1:length(fn)
        temp = jump_Ladder(FS_snpn1(i,1)).(fn{j})(:,FS_snpn1(i,2)) ;
        stack_jl1(ind1_1(i):ind2_1(i), j) = temp;
    end 
    clear temp
    for j = 1:length(fnR)
        temp = jump_Rat(FS_snpn1(i,1)).(fnR{j})(:,FS_snpn1(i,2)) ;
        stack_jr1(ind1_1(i):ind2_1(i), j) = temp;
    end 
    clear temp
    for j = 1:length(fn)
        temp = jump_Cat(FS_snpn1(i,1)).(fn{j})(:,FS_snpn1(i,2)) ;
        stack_jc1(ind1_1(i):ind2_1(i), j) = temp;
    end  
    
    stack_Npre1(ind1_1(i):ind2_1(i), 1) = seg_CN_pre{FS_snpn1(i,1)}(:,FS_snpn1(i,2));
    stack_Npre1(ind1_1(i):ind2_1(i), 2) = seg_CS_pre{FS_snpn1(i,1)}(:,FS_snpn1(i,2));
%     
%     clear temp
%     for j = 1:20
%         temp = broken_Trigger{FS_snpn1(i,1)}{FS_snpn1(i,2),j} ;
%         stack_broke1(ind1_1(i):ind2_1(i), j) = temp;
%     end 
    
end
stack_cc1 = [(stack_bc1(:,1)+stack_bc1(:,4))/2 , ...
    (stack_bc1(:,2)+stack_bc1(:,5))/2, (stack_bc1(:,3)+stack_bc1(:,6))/2];
stack_jr1(isnan(stack_jr1)) = 0;
stack_jr1(isinf(stack_jr1)) = 0;
% cross correlation between each fields, bdx, bdy, ... cn, cs, but only for
% those located with x>=1cm
jl_cor1 = corrcov(cov(normalize(stack_jl1(stack_bc1(:,1)>=0.01,:),'range')));
jc_cor1 = corrcov(cov(normalize(stack_jc1(stack_bc1(:,1)>=0.01,:),'range')));
jr_cor1 = corrcov(cov(normalize(stack_jr1(stack_bc1(:,1)>=0.01,:),'range')));                                                                                                                                                                                                                                                                                                                                                                                         jr_cor1 = corrcov(cov(normalize(stack_jr1,'range')));
%
[r, LB, UB, F, df1, df2, p] = ICC(stack_jr1(stack_bc1(:,1)>=0.01,:)', '1-1');
%
ind1_2=zeros(size(FS_snpn2,1),1); ind2_2=zeros(size(FS_snpn2,1),1); leng_2=zeros(size(FS_snpn2,1),1);
for i = 1:size(FS_snpn2)
    sn = FS_snpn2(i,1);
    leng_2(i) = size(Bcoord{sn},1);
    if i == 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
        ind1_2(i) = 1+0; ind2_2(i) = ind1_2(i)+leng_2(i)-1;
    else
        ind1_2(i) = 1+ind2_2(i-1); ind2_2(i) = ind1_2(i)+leng_2(i)-1;
    end     
end
mem_alloc_rows_2 = ind2_2(end);
%
stack_jl2 = zeros(mem_alloc_rows_2, length(fn));
stack_jc2 = zeros(mem_alloc_rows_2, length(fn)); 
stack_jr2 = zeros(mem_alloc_rows_2, length(fnR));
stack_bc2 = zeros(mem_alloc_rows_2, 6);
stack_cc2 = zeros(mem_alloc_rows_2, 3);
for i = 1:size(FS_snpn2,1)
    for j = 1:length(fn)
        temp = jump_Ladder(FS_snpn2(i,1)).(fn{j})(:,FS_snpn2(i,2)) ;
        stack_jl2(ind1_2(i):ind2_2(i), j) = temp;
    end
    clear temp
    temp = Bcoord{FS_snpn2(i,1)};
    stack_bc2(ind1_2(i):ind2_2(i), :) = temp;
% fix this! change to FS_snpn2    
    clear temp
    temp = Bcoord{FS_snpn2(i,1)};
    stack_bc2(ind1_2(i):ind2_2(i), :) = temp;
    clear temp
    for j = 1:length(fn)
        temp = jump_Ladder(FS_snpn2(i,1)).(fn{j})(:,FS_snpn2(i,2)) ;
        stack_jl2(ind1_2(i):ind2_2(i), j) = temp;
    end 
    clear temp
    for j = 1:length(fnR)
        temp = jump_Rat(FS_snpn2(i,1)).(fnR{j})(:,FS_snpn2(i,2)) ;
        stack_jr2(ind1_2(i):ind2_2(i), j) = temp;
    end 
    clear temp
    for j = 1:length(fn)
        temp = jump_Cat(FS_snpn2(i,1)).(fn{j})(:,FS_snpn2(i,2)) ;
        stack_jc2(ind1_2(i):ind2_2(i), j) = temp;
    end 
    %
end
stack_cc2 = [(stack_bc2(:,1)+stack_bc2(:,4))/2 , ...
    (stack_bc2(:,2)+stack_bc2(:,5))/2, (stack_bc2(:,3)+stack_bc2(:,6))/2];
%% keyword stacking - Part II
%ind1{3}

%temp = 0:0.25e-5:3.5e-5;
%temp = 0:0.1e-5:2.5e-5;
%temp = 7.75e5:-0.25e5:0e5; % For S200 series of simulations
%temp = 6e5:-0.5e5:0e5; % For S400 series of simulations
FS_split = 6e5:-0.25e5:0e5;
temp = FS_split;
FS_snpn = cell(length(temp)+4-1, 1);
nfssnpn = length(FS_snpn);
bnds_fssnpn = temp;
FS_snpn{1} = FS_snpn1;
FS_snpn{2} = FS_snpn2;
FS_snpn{3} = FS_snpn3;
FS_snpn{4} = FS_snpn0;%during window_stickslip
stack_fskdfp = zeros(length(temp)-1,1);
for i = 2:length(temp)
%    FStestgroupidx =  find( FStotal(:,2)<temp(i) & FStotal(:,2)>temp(i-1));%  & ...
%         FStotal(:,1)<3e5); %FStotal(:,3)<mv_W/4 & FStotal(:,3)>-mv_W/4 & ...
    %FStestgroupidx =  find( FStotal(:,2)<temp(i) ); 
    FStestgroupidx = find( FStotal(:,1)>temp(i) & FStotal(:,1)<temp(i-1));
    stack_fskdfp(i-1) = mean(FStotal(FStestgroupidx, 2));
    FS_snpn{i+3} = FStotal(FStestgroupidx, 6:7);
end

stack_jl = cell(nfssnpn,1);
stack_jc = cell(nfssnpn,1);
stack_jr = cell(nfssnpn,1);
stack_jpre = cell(nfssnpn,1);
stack_jpost = cell(nfssnpn,1);
stack_bc = cell(nfssnpn,1);
stack_cc = cell(nfssnpn,1);
stack_Npre = cell(nfssnpn,1);


for ii = 1:nfssnpn
temp = zeros(size(FS_snpn{ii},1),3);
ind_{ii} = temp;
for i = 1:size(FS_snpn{ii})
    sn = FS_snpn{ii}(i,1);
    ind_{ii}(i,3) = size(Bcoord{sn},1);
    if i == 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
        ind_{ii}(i,1) = 1+0; ind_{ii}(i,2) = ind_{ii}(i,1)+ind_{ii}(i,3)-1;
    else
        ind_{ii}(i,1) = 1+ind_{ii}(i-1,2); ind_{ii}(i,2) = ind_{ii}(i,1)+ind_{ii}(i,3)-1;
    end     
end
mem_alloc_rows(ii) = ind_{ii}(end,2);
stack_jl3 = zeros(mem_alloc_rows(ii), length(fn));
stack_jc3 = zeros(mem_alloc_rows(ii), length(fn)); 
stack_jr3 = zeros(mem_alloc_rows(ii), length(fnR));
stack_bc3 = zeros(mem_alloc_rows(ii), 6);
stack_cc3 = zeros(mem_alloc_rows(ii), 3);
stack_Npre3 = zeros(mem_alloc_rows(ii), 3);
for i = 1:size(FS_snpn{ii},1)
    for j = 1:length(fn)
        temp = jump_Ladder(FS_snpn{ii}(i,1)).(fn{j})(:,FS_snpn{ii}(i,2)) ;
        stack_jl3(ind_{ii}(i,1):ind_{ii}(i,2), j) = temp;
    end
    clear temp
    temp = Bcoord{FS_snpn{ii}(i,1)};
    stack_bc3(ind_{ii}(i,1):ind_{ii}(i,2), :) = temp;   
%     clear temp
%     temp = Bcoord{FS_snpn2(i,1)};
%     stack_bc2(ind1_2(i):ind2_2(i), :) = temp;
    clear temp
    for j = 1:length(fn)
        temp = jump_Ladder(FS_snpn{ii}(i,1)).(fn{j})(:,FS_snpn{ii}(i,2)) ;
        stack_jl3(ind_{ii}(i,1):ind_{ii}(i,2), j) = temp;
    end 
    clear temp
    for j = 1:length(fnR)
        temp = jump_Rat(FS_snpn{ii}(i,1)).(fnR{j})(:,FS_snpn{ii}(i,2)) ;
        stack_jr3(ind_{ii}(i,1):ind_{ii}(i,2), j) = temp;
    end 
    clear temp
    for j = 1:length(fn)
        temp = jump_Cat(FS_snpn{ii}(i,1)).(fn{j})(:,FS_snpn{ii}(i,2)) ;
        stack_jc3(ind_{ii}(i,1):ind_{ii}(i,2), j) = temp;
    end 
    stack_Npre3(ind_{ii}(i,1):ind_{ii}(i,2), 1) = seg_CN_pre{FS_snpn{ii}(i,1)}(:,FS_snpn{ii}(i,2));
    stack_Npre3(ind_{ii}(i,1):ind_{ii}(i,2), 2) = seg_CS_pre{FS_snpn{ii}(i,1)}(:,FS_snpn{ii}(i,2));
    
    %
end
stack_Npre3(:,3) = stack_Npre3(:,1)/(2e-5)*tand(30)+55e6;%estimated shear strength
stack_cc3 = [(stack_bc3(:,1)+stack_bc3(:,4))/2 , ...
    (stack_bc3(:,2)+stack_bc3(:,5))/2, (stack_bc3(:,3)+stack_bc3(:,6))/2];
stack_jl{ii} = stack_jl3;
stack_jc{ii} = stack_jc3; 
stack_jr{ii} = stack_jr3;
stack_bc{ii} = stack_bc3;
stack_cc{ii} = stack_cc3;
stack_Npre{ii} = stack_Npre3; %1st column normal force, 2nd column shear force, 3rd column shear strength
stack_jpre{ii}= stack_jl{ii}./stack_jr{ii}(:,1:11);
stack_jpost{ii}= stack_jl{ii}./stack_jr{ii}(:,12:22);

stack_jpre{ii}(isnan(stack_jpre{ii})) = 0;
stack_jpost{ii}(isnan(stack_jpost{ii})) = 0;
end


%% keyword hist
edges = linspace(0,3e8,100);%0:0.2e8:3e8;
hist_oa = zeros(nfssnpn, 1);
hist_oa2 = cell(nfssnpn,1);
for i = 1:nfssnpn
    figure1 = figure(f);f=f+1;
    h1 = histogram(nonzeros(stack_Npre{i}(:,2)/(2e-5)), 'BinEdges',edges,...
        'Normalization', 'pdf'); %histogram of shear stress
    hold on;
    %h2 = histogram(stack_Npre{i}(stack_Npre{i}(:,3)>55e6,3), 'BinEdges',edges,...
    %    'Normalization', 'pdf');hold off; %histogram of shear strength that are greater than cohesion
    h2 = histogram(nonzeros(stack_Npre{i}(:,3)), 'BinEdges',edges,...
        'Normalization', 'pdf');hold off; %histogram of shear strength
    ylabel('Probability density');
    legend('shear stress','shear strength');
    set(gca,'FontSize', 12,'FontName', 'Times');
    figurename = strcat('400stack_Npre',num2str(i),'.png'); saveas(figure1,figurename);
    
    temp1 = h1.BinEdges; temp2 = h2.BinEdges;
    imin = find(temp1==min(temp2));
%     hist_oa(i) = sum(min( ...
%         h1.Values(imin:length(temp1)-1),h2.Values(1:length(temp1)-imin)...
%         ))*h1.BinWidth; %histogram overlapping area
    
%     hist_oa(i) = sum( ...
%         h1.Values(imin:length(temp1)-1)...
%         )*h1.BinWidth; %histogram overlapping area
    temp3 = h1.Values; temp4 = h2.Values;
    hist_oa(i) = sum( min(temp3(temp4~=0&temp3~=0), ...
        temp4(temp4~=0&temp3~=0)) )*h1.BinWidth; 
%     hist_oa2{i} = (temp4(temp4~=0&temp3~=0)-temp3(temp4~=0&temp3~=0))...
%         *h1.BinWidth;
    hist_oa2{i} = (temp4(edges(1:end-1)>55e6)-temp3(edges(1:end-1)>55e6))...
        *h1.BinWidth;

%     hist_dist(i) = sum(...
%         temp3(temp4~=0&temp3~=0).*...
%         log(temp3(temp4~=0&temp3~=0)./temp4(temp4~=0&temp3~=0))...
%         .*h1.BinWidth...
%         ); %Kullback–Leibler divergence, does not work
%     
    %sum(h2.Values)*h2.BinWidth - sum(h1.Values(imin:length(temp1)-1))*h1.BinWidth
    
    %clear temp1; clear temp2;
end
%
close all
fig=figure(f);f=f+1;scatter(-bnds_fssnpn(1:end-1),hist_oa(5:end), 'k');
%figure(f);f=f+1;scatter(bnds_fssnpn,[hist_oa(5:end);hist_oa(4)]);
%xlabel('time to mainshock (sec)');
xlabel('Time relative to the onset of fault slip [steps]')
%xlabel('KDF peak heights');
%ylabel('overlapping area under normalized histograms');
ylabel({'Intersection of', 'shear stress and strength pdfs'});
set(gca,'FontSize', 12, 'FontName', 'Times');
%exportgraphics(fig,'evol of overlapping area.jpg','Resolution',300);



%% (keyword hist) spatial distribution of jump ladder - difference of areas under histogram of shear and normal forces
% 3D histogram, each slice along y direction is a histogram
%temp = 6e5:-0.5e5:0e5;
meow = zeros(length(FS_split)+2-5,1);
nedges = 50;%100;
edges = linspace(-8e-5,8e-5,nedges);%linspace(-8e-5,8e-5,50);
histn = zeros(nedges-1, 2);
%fff = [15, 30];
fff = 5:length(FS_split)+2; 
diffmodelever=zeros(length(FS_split)+2,2,3); diffavg=zeros(length(FS_split)+2,2,3);diffmode=zeros(length(FS_split)+2,2,3);
for ff = 1:length(fff) %fi controls which time segment
    %figure1 = figure(f);f=f+1;
    for k = [10,11]%[1,3,5,10,11] %k controls which fn, [1,3,5,10,11], 10 for CN, 11 for CS
        temp1 = stack_jpre{fff(ff)}(:,k);
        temp2 = stack_jpost{fff(ff)}(:,k);
        temp3 = temp2 - temp1;
        %temp3 = rmoutliers(temp3);
        %h1 = histogram(temp3, 'Normalization','probability', 'BinEdges', edges);title('distrbtn. of force accel. Orange for shear & blue for normal')
        h1 = histcounts(temp3, 'Normalization','probability', 'BinEdges', edges);
        hold on
        histn(:,k-9) = h1;
    end
    hold off; 
    %saveas(figure1,strcat(fgroup,'dist_f_accel',num2str(ff),'.png'));
    meow(ff,1) =sum( histn(:,2) - histn(:,1))*(edges(2)-edges(1));
end

x = (FS_split(1:end-1) + FS_split(2:end))/2;


fig = figure(f);f=f+1;
scatter(x(11:end-1)*-1,meow(11:end), 'k'); 
ylabel({'Evolution of overall difference'; 'of ROC between';'shear and normal forces [N/step]'});
xlabel('Time relative to the onset of fault slip [steps]');
set(gca,'FontSize', 12, 'FontName', 'Times');
exportgraphics(fig,'evol of ROC.jpg','Resolution',300);

