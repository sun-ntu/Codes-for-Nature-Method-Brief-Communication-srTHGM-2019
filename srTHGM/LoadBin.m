function A = LoadBin(Folder,Ch,y,r,ZL,z)

    A = [];
    
    subFolder = [Folder '/Channel_' num2str(Ch)];
    allfile = dir(fullfile(subFolder,'*Channel*'));
    allfile_new=dir(fullfile(subFolder,'*CHAN_*'));
    if numel(allfile_new)==0
        A = ZL-Bin2img([subFolder '/' num2str(Ch) '_Channel_' num2str(z)],y,r);
    elseif numel(allfile)==0
        temp = Bin2img_stk([subFolder '/' 'CHAN_' num2str(Ch)],y,r);
        A=temp(:,:,z);
        A=ZL-A;
    else
        error('No supported file format')
    end
    
end