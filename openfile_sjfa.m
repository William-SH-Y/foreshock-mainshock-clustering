function [num_record, b_output] = openfile_sjfa(fname,num_sensor,num_rec)
    %filename = fopen('bvxS51.txt');
    filename = fopen(fname);
    tempstring = '%f';

    for i = 1:num_sensor
        tempstring = strcat(tempstring, ' %f');
    end
    C = textscan(filename, tempstring);
    bvtemp = zeros(num_rec, num_sensor);
    for i = 1:num_sensor
        bvtemp(:, i) = cell2mat(C(1,i+1));
    end
    num_record = nnz(bvtemp(:, 1));
    b_output = zeros(num_record, num_sensor);
    for i = 1:num_sensor
        b_output(:, i) = bvtemp(1:num_record, i);
    end
    clear bvtemp
    clear tempstring
end