function W = Solve_W2_new(E, P,mu_B)
% solve the following problem:
% mu_B/2*||W-E||^2_2+k*||P*(卷积)W||^p_p
% 转化为：mu_B/2*||W-E||^2_2+ita*k/2*||P*(卷积)W-d||^2_2 和 ita/2*||d-d0||^2_2+||d||^p_p 求解，
% 其中 d0=P*(卷积)E
q        =  0.4;  
tau      =  1;
%  a        =  1;            %a = ita*tau/mu_B   
ita      =  0.02;         %ita控制平滑程度，ita越小平滑程度越大（试过0.02和0.2）
rho      =  1.05;         %rho = 2*sqrt(2)

% 在此处调节迭代次数：判断ASR的效果！
% 对应期刊上的图3： k = 5, 10, 15, 20
itermax  =  20;  % 50 

all_psnr =  zeros(1,itermax);
K        =  size(P,1);
[m,n]    =  size(E);
% 初始化
W  = E;
% FP = zeros(m,n,K,'single');
FP = zeros(m,n,K);
for i=1:K
    FP(:,:,i) = psf2otf(reshape(P(i,:),[8,8]),[m,n]);
    %FP(:,:,i) = psf2otf(reshape(P(i,:),[3,3]),[m,n]);
end
FW = repmat(fft2(W),[1,1,K]);
Fd = FP.*FW;
FE = fft2(E);
SumFPTP = sum(FP.*conj(FP),3);

% while ita<1e+6     %ita<2^8
    for iter = 1:itermax
        a = ita*tau/mu_B;
        FW    =  (FE+a*sum(conj(FP).*Fd,3))./(ones(size(E))+a*SumFPTP);   
        W     =  real(ifft2(FW));
        
        d0    =  real(ifft2(FP.*repmat(FW,[1,1,K])));
        d     =  GST( d0, 1/ita, q);
        Fd    =  fft2(d);
        
%         a     =  rho*a;
        ita   =  rho*ita;
%         figure(1);imshow(uint8(W))
%         all_psnr(iter)  =  csnr(W,E,0,0);
%       if iter>1
%         if(all_psnr(iter) - all_psnr(iter-1)<0)
%             break;
%         end
%       end
    end
% end

