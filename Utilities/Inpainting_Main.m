function [corrupt_image , x_inpaint_rgb_all , All_PSNR]=Inpainting_Main(ori_gname,patch,p_miss,delta,mu_A,mu_B,mu_C,p)
        
        x_rgb                  =    imread(ori_gname); 
        
        x_yuv                  =    rgb2ycbcr(x_rgb);
        
        x                      =    double(x_yuv(:,:,1)); 
        
        x_org                  =    x;
        
        x_inpaint_re_all          =    zeros(size(x_yuv));
        x_inpaint_re_all(:,:,2)    =    x_yuv(:,:,2); 
        x_inpaint_re_all(:,:,3)    =    x_yuv(:,:,3); 
        
        ratio                  =    p_miss; 


        rand('seed',0);
        O = double(rand(size(x)) > (1-ratio));

        randn('seed',0);
        x = x + delta*randn(size(x));
        y = double(x).* O;       % Observed Image
        
     
        para = [];
        
        if ~isfield(para,'mu_A')
            para.mu_A = mu_A;
        end
        
        if ~isfield(para,'mu_B')
            para.mu_B = mu_B;
        end

        if ~isfield(para,'mu_C')
            para.mu_C = mu_C;
        end
        
        if ~isfield(para,'org')
            para.org = x_org;
        end  
        
        if ~isfield(para,'IterNums')
            para.IterNums = 200;  
        end 
        
        if ~isfield(para,'Initial')
            para.Initial = Inter_Initial(y,~O);
%             para.Initial = y;
        end
        
        if ~isfield(para,'patch')
            para.patch = patch;  % patch d
        end
        
        if ~isfield(para,'step')
            para.step = 2;
        end       
        
        if ~isfield(para,'Similar_patch')
            para.Similar_patch = 60; % Similar patches c
        end
         
        if ~isfield(para,'Region')
            para.Region = 25;
        end        
        
        if ~isfield(para,'sigma')
            para.sigma = sqrt(2);
        end 
        
        if ~isfield(para,'e')
            para.e = 0.3;
        end         
        
%         if ~isfield(para,'nSig')
%             para.nSig = 30;
%         end  
               
     fprintf(ori_gname);
     fprintf('\n');

     [reconstructed_image_all, All_PSNR]  =  Inpainting_all(y,O,para,p); 


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     corrupt_image =  zeros(size(x_yuv));
     corrupt_image(:,:,1) = y;
     corrupt_image(:,:,2) = x_yuv(:,:,2);
     corrupt_image(:,:,3) = x_yuv(:,:,3);
     corrupt_image = ycbcr2rgb(uint8(corrupt_image));
     
     x_inpaint_re_all(:,:,1) = reconstructed_image_all;
     x_inpaint_rgb_all = ycbcr2rgb(uint8(x_inpaint_re_all));
         
     
     fprintf('Ending Algorithm for Image Inpainting\n'); 
     fprintf('................................................\n');       
       
end

