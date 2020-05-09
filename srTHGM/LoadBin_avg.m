function A = LoadBin_avg(Folder,Ch,y,r,ZL)

    A = 0;
    
    subFolder = [Folder '/Channel_' num2str(Ch)];
    allfile = dir(fullfile(subFolder,'*Channel*'));
    allfile_new=dir(fullfile(subFolder,'*CHAN_*'));
    if numel(allfile_new)==0
%         for i2 = 1:numel(allfile)  % for conveience
         for i2 = 1:3  % for conveience
            temp = Bin2img([subFolder '/' allfile(i2).name],y,r);
            A = A+temp;
        end
%         A = ZL-A/numel(allfile);
        A = ZL-A/3;
        %     A = ZL-Bin2img([subFolder '/' allfile(1).name],y,r);
    elseif numel(allfile)==0
        temp = Bin2img_stk([subFolder '/' 'CHAN_' num2str(Ch)],y,r);
         temp = temp(:,:,4:6);
        temp_avg=mean(temp,3);
        A=ZL-temp_avg;
    else
    end
end