function A = Bin2img_stk(filename,y,r)

fileID = fopen(filename, 'r', 'ieee-le');
if fileID == -1, error('Cannot open file: %s', filename); end
format = 'uint16';
Data = fread(fileID, Inf, format);
fclose(fileID);
z = length(Data)/3024/y;
% size(Data)
A = reshape(Data,3024,y,[]);

if r==1
    A = permute(A,[2 1 3]);
end