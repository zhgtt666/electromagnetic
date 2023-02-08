function [] = savegif(it,filename,sliceInterval,snpshotDelayTime,isinitial)
% Description: make gif snapshot file
% Version: 1.0
% INPUT
% it: the serial number of snapshot
% filename: file name
% sliceInterval: every 'sliceInterval' slice steps save a frame
% snapshotDelayTime: two frame interval time in gif
%isinitial: =0£ºadd frame into gif file£¬=others£ºinitialize the gif file 
% Note: when first time call this function, it must be equal to 1
% ------------------
% Autor: Chen QuYang
% Date: 2022-12-31
% LastEditors: Chen QuYang
% LastEditTime: 2023-01-13
% Copyright (c) 2023 WaveTomo. All rights reserved. 
frame = getframe(gcf);
imind = frame2im(frame); 
[imind, cmap] = rgb2ind(imind, 256); 
if isinitial ~=0 
    imwrite(imind,cmap, filename,'gif', 'Loopcount',inf,'DelayTime',snpshotDelayTime);
else 
    if(mod(it,sliceInterval)==0)
        imwrite(imind,cmap, filename,'gif','WriteMode','append','DelayTime',snpshotDelayTime);
    end
end
end

