function A = Bin2img(filename,y,r)

fileID = fopen(filename, 'r', 'ieee-le');
if fileID == -1, error('Cannot open file: %s', filename); end
format = 'uint16';
Data = fread(fileID, Inf, format);
fclose(fileID);

A = reshape(Data,[],y);

if r==1
    A = A';
end