function  [sTimb,sTim] = blurimN(sT,rM,Mroi,u,cmap2)

r = u+2;
r2 = (r+1)/2;
u2 = (u+1)/2;
ag = zeros(r+2,r+2);

for ix=1:r
    for iy = 1:r
        
        t = (ix-r2)^2+(iy-r2)^2;

        if(t<=(u2-0.5)^2)
            ag(ix,iy) = exp(-t/r);
        end
    end
end
ag = ag/max(max(ag));

sT = sT.*(sT>0.01)+0.01*(sT<=0.01).*rM;

im = ind2rgb(uint16(sT*2^16),cmap2);
%%
% t4 = conv2(rM,ae,'same');

t4c = conv2(rM,nthroot(ag,(40/u)),'same');
t4c = t4c + 100*double(t4c==0);

%%

t6 = im;
t6(:,:,1) = conv2(t6(:,:,1),ag,'same')./t4c;
t6(:,:,2) = conv2(t6(:,:,2),ag,'same')./t4c;
t6(:,:,3) = conv2(t6(:,:,3),ag,'same')./t4c;

% figure(n+4)
% imshow(t6)
%%

kk = (t4c==100);

sTimb = zeros(size(t6));
sTimb(:,:,1) = t6(:,:,1).*Mroi+0.25*Mroi.*kk;
sTimb(:,:,2) = t6(:,:,2).*Mroi+0.25*Mroi.*kk;
sTimb(:,:,3) = t6(:,:,3).*Mroi+0.25*Mroi.*kk;

% figure(n+5)
% imshow(t7)
%%

kk2 = (sT==0);

sTim = zeros(size(t6));
sTim(:,:,1) = im(:,:,1)+0.25*Mroi.*kk2;
sTim(:,:,2) = im(:,:,2)+0.25*Mroi.*kk2;
sTim(:,:,3) = im(:,:,3)+0.25*Mroi.*kk2;

