function data_converted_fs = convert_ivs_to_other_fs(data,curr_fs, new_fs)
data_seconds = data/curr_fs; %convert to seconds
data_converted_fs = round(data_seconds*new_fs); %multiply to new fs
end
%data*new_fs/currfs