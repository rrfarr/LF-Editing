function [lf_name, ds_name, psnr, ssim] = csv_read_SR_quality_file(fid)

k = 1;
while ~feof(fid)
    % Read a line from the csv file
    line = fgetl(fid);
    % Derive the index of comma
    idx_comma = strfind(line,',');
    lf_name{k,1} = line(1:idx_comma(1)-1);
    ds_name{k,1} = line(idx_comma(1)+1:idx_comma(2)-1);
    psnr(k,1)    = str2double(line(idx_comma(2)+1:idx_comma(3)-1));
    ssim(k,1)    = str2double(line(idx_comma(3)+1:end));
    k = k + 1;
end




