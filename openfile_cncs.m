function [b_output] = openfile_cncs(fname,num_sensor, num_rec, num_record)

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

    b_output = zeros(num_record, num_sensor);
    for i = 1:num_sensor
        b_output(:, i) = bvtemp(1:num_record, i);
    end
    clear bvtemp
    clear tempstring
end