clc;
close all;
clear all;
%% Read the Image
Image = im2double(imread('BME7112_Data_File_12.tif'));
% used resize to perform fast operation in the loop
% Image = imresize(Image,[1000 1000]);
% Image= uint8(Image / 256);
Image1 = adapthisteq(Image,'ClipLimit',0.1,'NumTiles',[4,4]);% histograam equlization is used to enhance the contrast in the image
% SO while cropping it is easily visible.
figure('Name','Raw Image')
imshow(Image1)
title('Contrast enhanced raw image')
help = helpdlg('Crop the intrested metacarpals');% this help dialoge box will remind user to choose the ROI, here user has to 
%click ok after clicking ok, the cropping will initiative
uiwait(help)
[Image1,position] = imcrop(Image1);% cropping the Image . This is the user input to select the intrested metacarpals
figure('Name','Cropped Image')
imshow(Image1)
title('Contrast enhanced cropped raw image')
% hRoi = imfreehand();% make freehand ROI From the rect ROI to analyse the focused part from the  Rect ROI
% Position = getPosition(hRoi);% finding the X and Y vectors of Freehand ROI
% BW = createMask(hRoi);% MAking a mask 
% A(BW == 0) = 0;
figure('Name','Cropped contrast enhanced Image')
Image = imcrop(Image,position);% need to use histogram equilization on the cropped image, SO we can enhance much
% better enhanced image in the region of intrest
% Image = Image.*BW;
Image = adapthisteq(Image,'ClipLimit',0.1,'NumTiles',[4,4]);% straching the contrast in the ROI
% Image = imsharpen(Image,'Radius',4,'Amount',4);
imshow(Image)
title('after cropped contrast enhanced image')

%% Define the kernel size an padding the image 
Kernel_size = 9;% define the kernel size, by the reference paper, Kernel size has to 9.
Image = padarray(Image,[Kernel_size Kernel_size],0,'both');% padding the image with 0 for both side same as kernel size, to do efficient 
%sliding window operation on image
Image_core = Image;% image core variable is used to ovaelay the segmented mask on the Image_core variable
[x,y] = size(Image);% finding the size of the image to perform sliding window operation
% Imagex1 = zeros(x,y);

%% sliding window operation to perform the praposed algorhytm
start = (round(Kernel_size/2));% for odd kernel my start point is middle pixel of the kernel so this is the adjustmnt to find middele pixel
as = Kernel_size - start+1;% adjustment to find the central pixel
filtered_ImageA = zeros(x,y);  % the updating pixel in the loop will replaced in this matrix so original matrix values are not chnaged.
% this loop suits the odd number kernel size not even.
    for i = start:x-start+1
        for j = start:y-start+1
            sum = 0;
            for k = 1:Kernel_size
                for l = 1:Kernel_size
                    sum(k,l) = Image(k+i-as,l+j-as);
                end
            end
              Kernel_Entropy = (wentropy(sum,'shannon'));% finding the shannon entropy
              Kernel_Std = std(sum(:));% finding standard deviation 
              Kernel_Mean = mean(sum(:));% finding mean
              filtered_ImageA(k+i-Kernel_size,j+l-Kernel_size) = Kernel_Entropy*Kernel_Std/Kernel_Mean;% parpoped algorhytm will chnage he central pixel value with the e value.
%             T = graythresh(filtered_ImageA(i:i+k-1,j+l-1));try to
%             threshold it in the loop but no working
% %             BW = imbinarize(filtered_ImageA(i:i+k-1,j:j+l-1),T);
%             Filtered_ImageB(i:i+k-1,j:j+l-1) = e;
%              BW = imbinarize(filtered_ImageA(k+i-Kernel_size,j+l-Kernel_size),T);
%              Filtered_ImageB(k+i-Kernel_size,j+l-Kernel_size) = e;

        end 
    end
% filtered_ImageA(isnan(filtered_ImageA))=0;
%% Normalization of the filtered image
filtered_ImageA_Max  = max(filtered_ImageA(:));% finding the max value from the filtered_image
filtered_ImageA = filtered_ImageA./filtered_ImageA_Max;% normlaizing the filtered_Image between [0 , 1]

%% ploting the filtered image
figure('Name','entropy image')
imshow(filtered_ImageA)
title('Normlaized image after praposed algorhytm')
filtered_ImageA = imcomplement(filtered_ImageA);% comlement the image to see sorrounding tissue brighter than bone
figure('Name','Complement of entropy image')
imshow(filtered_ImageA)
title('Complement of normalized image')
%   filtered_ImageA = adapthisteq(filtered_ImageA,'ClipLimit',0.1,'NumTiles',[4,4]);
%   filtered_ImageA = imsharpen(filtered_ImageA,'Radius',2,'Amount',2); 
% try to sharp and enhanced image after filtering but that enhanced my
% background noise too much so I have used unsharp mask
%% unsharp mask to remove the background
filtered_ImageA_Smooth = imgaussfilt(filtered_ImageA,10);% used gaussian filter to smooth out the praposed algorhytm.
filtered_ImageA_edge_enhanced = filtered_ImageA_Smooth - filtered_ImageA; % 
figure('Name','edge enhance image with minimal noise')
imshow(filtered_ImageA_edge_enhanced)
title('Image after removing background and enhance the edge')
%% thresholding the filtered image

T = graythresh(filtered_ImageA_edge_enhanced);% finding threshold value according to the otsu's threhsolding method.
BW = imbinarize(filtered_ImageA_edge_enhanced,T);% binarize the imahe according to threshold
figure('Name','threshold image')
imshow(BW)
title('BW mask after thresholding')

%% performing morphological operations

bw1 = bwmorph(BW,'skel',Inf);% perform skeltization because my Bw mask has thick edges so this will make it more clear 
%as recommended(1 pixel wide)
figure('Name','Skel')
imshow(bw1)
title('BW mask after skeletization')
% BW1 = bwmorph(bw1,'clean',Inf);
% figure()
% imshow(BW1)
% title('after clean')
% 
% BW2 = bwmorph(BW1,'close',Inf);
% figure()
% imshow(BW2)
% title('after close')
% 
% BW3 = bwmorph(bw1,'remove',Inf);
% figure()
% imshow(BW3)
% title('after remove')
% 
% BW4 = bwmorph(BW3,'thin',Inf);
% figure()
% imshow(BW4)
% title('after thin')

% BW5 = bwmorph(BW4,'skel',Inf);
% figure()
% imshow(BW6)
% title('Skel')
% title('open image')
BW6 = bwareaopen(bw1,20);% perform this operation to remove unwanted edges from the background as well as in the bone.
figure('Name','bwareaopen')
imshow(BW6)
title('BW mask after skel & bwareaopen')

% BW7 = bwmorph(BW6,'skel',Inf);
% figure()
% imshow(BW7)
WE = imoverlay(Image_core,BW6,'green');
figure('Name','Overlayed Image')
imshow(WE)
title('Overlayed Image')
%%
% h = ones(Kernel_size,Kernel_size) / Kernel_size*Kernel_size;
% filtered_ImageA_avg = imfilter(filtered_Imageam,h);
% filtered_ImageA_avg_max = max(filtered_ImageA_avg(:));
% filtered_ImageA_avg_norm = filtered_ImageA_avg./filtered_ImageA_avg_max;
% figure()
% imshow(filtered_ImageA_avg_norm)
% title('Image After Averaging')
% final_ImageA = (Image - filtered_ImageA_avg_norm);
% % final_ImageA = abs(final_ImageA);
% final_ImageA(final_ImageA<0) = 0;
% final_ImageA(isnan(final_ImageA))=0;
% % final_ImageA = imcomplement(final_ImageA);
% % final_ImageA = adapthisteq(final_ImageA,'ClipLimit',0.1,'NumTiles',[4,4]);
% % final_ImageA = imsharpen(final_ImageA,'Radius',1,'Amount',1);
% T = graythresh(final_ImageA);
% % final_ImageA = edge(final_ImageA,'canny',T);
% % [counts,x] = imhist(final_ImageA,255);
% % T = otsuthresh(counts);
% % T =adaptthresh(final_ImageA, 0.2);
% BW1 =  imbinarize(final_ImageA,T);
% % BW1 = ~BW1;
% % BW2 = edge(BW1,'sobel');
% % figure()
% % imshow(BW2)
% % title('edge detection segment after averaging')
% % figure()
% % imshow(BW1)
% % title('After subtracting from Averaging image')
% % bw1 = bwmorph(BW1,'skel',Inf);
% % BWq = bwmorph(bw1,'clean',Inf);
% % BW2 = bwmorph(BWq,'close',Inf);
% % BW3 = bwmorph(BW2,'remove',Inf);
% % BW4 = bwmorph(BW3,'thin',Inf);
% % WE = imoverlay(Image_core,BW4,'green');
% % figure()
% % imshow(WE)
% % figure()
% % imshow(final_ImageA)
% % title('final Image')
% bw1 = bwmorph(BW1,'skel',Inf);
% figure()
% imshow(bw1)
% title('after skel')
% 
% BW6 = bwareaopen(bw1,20);
% figure()
% imshow(BW6)
% WE = imoverlay(Image_core,BW6,'green');
% figure()
% imshow(WE)
