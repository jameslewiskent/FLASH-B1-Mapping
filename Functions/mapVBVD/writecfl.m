function writecfl(filenameBase,data)
%WRITECFL  Write complex data to file.
%   Writes reconstruction data to filenameBase.cfl (complex float) and its
%   dimensions to filenameBase.hdr.
%
%   Written to edit data with the Berkeley Advanced Reconstruction Toolbox (BART).
%
% Copyright 2013. Joseph Y Cheng.
% Copyright 2016. CBClab, Maastricht University.
% 2012 Joseph Y Cheng (jycheng@mrsrl.stanford.edu).
% 2016 Tim Loderhose (t.loderhose@student.maastrichtuniversity.nl).

    dims = size(data);
    writeReconHeader(filenameBase,dims);

    filename = strcat(filenameBase,'.cfl');
    fid = fopen(filename,'w');
    
    data = data(:);
    
    % ATH June 2017 for large data split it
    part_len = 2e9;
    file_len = length(data); %*2 real imaginary
    file_written = 0;
    if(file_len>part_len)
        while((file_len-file_written)>0)
            if((file_len-file_written)<part_len)
                part_len = file_len-file_written;
            end
            file_written = file_written + fwrite(fid,[real(data(file_written+1:part_len+file_written))'; imag(data(file_written+1:part_len+file_written))'],'float32')/2;
            
        end
    else
        fwrite(fid,[real(data)'; imag(data)'],'float32');
    end
    fclose(fid);
end

function writeReconHeader(filenameBase,dims)
    filename = strcat(filenameBase,'.hdr');
    fid = fopen(filename,'w');
    fprintf(fid,'# Dimensions\n');
    for N=1:length(dims)
        fprintf(fid,'%d ',dims(N));
    end
    if length(dims) < 5
        for N=1:(5-length(dims))
            fprintf(fid,'1 ');
        end
    end
    fprintf(fid,'\n');
    
    fclose(fid);
end

