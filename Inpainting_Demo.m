% the proposed image inpainting method: CSRNS, which is submitted to SMCA
% a demo on image inpainting under noise level 30 with 50% missing pixel 

clc;
clear;
close all;
addpath('./Test_Images/');
addpath('./Utilities/');




for i = 1 : 7
    
    for ImgNo = i
        switch ImgNo
            case 1
                filename = 'Butterfly';
            case 2
                filename = 'Castle';
            case 3
                filename = 'Corn';
            case 4
                filename = 'Mickey';
            case 5
                filename = 'Girl';
            case 6
                filename = 'House';
            case 7
                filename = 'Barbara';

        end
        
        orgname = [filename '.png'];
        

        delta = 30;  % noise level 30
        p_miss   =  0.5; % 50% random pixels missing
        
      
        fprintf('ImgNo = %f\n',ImgNo)
        [corrupt_image , im_out_all , All_PSNR]  =  Inpainting_Main(orgname, 8, p_miss,delta,0.1,0.1,1.7,0.95);
        imshow([corrupt_image, im_out_all],[0,1]);

        
        
    end
end

