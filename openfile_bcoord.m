function [b_output] = openfile_bcoord(fname,num_sensor,num_ball_sets)
    filename = fopen(fname);
    tempstring = '%f';
    num_ball_sets = num_ball_sets * 3;
    for i = 1:num_ball_sets
        tempstring = strcat(tempstring, ' %f');
    end
    C = textscan(filename, tempstring);
    bvtemp = zeros(num_sensor, num_ball_sets);
    for i = 1:num_ball_sets
        bvtemp(:, i) = cell2mat(C(1,i));
    end
    num_record = nnz(bvtemp(:, 2));
    b_output = zeros(num_sensor, num_ball_sets);
    for i = 1:num_ball_sets
        b_output(:, i) = bvtemp(1:num_record, i);
    end
    clear bvtemp
    clear tempstring
end