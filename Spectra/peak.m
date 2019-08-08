function [cw,bd] = peak(t,z,x,y)

yn = y./max(y);
cw = x(round(mean(find(yn>0.9))));

temp = find(yn>0.494);
bd_1 = x(temp(1));
bd_2 = x(temp(end));
bd = abs(bd_2-bd_1);

% display(['   -> ' num2str(t) ' (' num2str(z) ')  FWHM: ' num2str(bd) '   /   CW: ' num2str(cw) ]);
