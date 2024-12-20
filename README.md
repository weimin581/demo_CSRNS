Title. Simultaneous Image Denoising and Completion  through Convolutional Sparse Representation and Nonlocal Self-similarity

link: https://www.sciencedirect.com/science/article/pii/S1077314224002972

Abstract. Low rank matrix approximation (LRMA) has been widely studied due to its capability of approximating original image from the degraded image. According to the characteristics of degraded images, image denoising and image completion have become research objects. Existing methods are usually designed for a single task. In this paper, focusing on the task of simultaneous image denoising and completion, we propose a weighted low rank sparse representation model and the corresponding efficient algorithm based on LRMA. The proposed method integrates convolutional analysis sparse representation (ASR) and nonlocal statistical modeling to maintain local smoothness and nonlocal self-similarity (NLSM) of natural images. More importantly, we explore the alternating direction method of multipliers (ADMM) to solve the above inverse problem efficiently due to the complexity of simultaneous image denoising and completion. We conduct experiments on image completion for partial random samples and mask removal with different noise levels. Extensive experiments on four datasets, i.e., Set12, Kodak, McMaster, and CBSD68, show that the proposed method prevents the transmission of noise while completing images and has achieved better quantitative results and human visual quality compared to 17 methods. The proposed method achieves (1.9\%, 1.8\%, 4.2\%, and 3.7\%) gains in average PSNR and (4.2\%, 2.9\%, 6.7\%, and 6.6\%) gains in average SSIM over the sub-optimal method across the four datasets, respectively. We also demonstrate that our method can handle the challenging scenarios well.


![ff1](https://github.com/user-attachments/assets/429a2a8e-650c-4808-9f1c-5f743988fde6)
Visual comparison of Butterfly mixed corrupted by additive white Gaussian noise and irregular scratches. (a) Ground Truth (PSNR/SSIM). (b) Input corrupted image (9.83/0.209). (c) WSNM (18.75/0.440). (d) LGSR (19.32/0.468). (e) SHI (24.52/0.835). (d) Proposed (24.96/0.848). Through comparison, our method achieves the performance closest to GT. The proposed method not only produces the most natural-looking effect but also effectively denoises the image and restores missing information simultaneously. Zoom-in is recommended to see the detailed comparisons.

![f0](https://github.com/user-attachments/assets/6e8bc86b-2c4e-4683-a7cd-de8058e855b1)
The flowchart illustrates our proposed method for simultaneous image denoising and completion. The degraded image is reconstructed using the proposed model by jointly incorporating the ASR and NLSM into existing WSNM to form the objective function Eq. (11), which is iteratively solved through the ADMM framework. Specifically, the objective function is divided into Z, W, and S sub-problems and solved separately, then the reconstruction result X is obtained by combining these three modules.
![bscd68](https://github.com/user-attachments/assets/8e554741-1ecb-4ac3-9c87-6017023e9d73)


![f9](https://github.com/user-attachments/assets/5b5589ee-9a86-484e-815b-1f12ae4e7391)

@article{yuan2024csrns, title={Simultaneous Image Denoising and Completion  through Convolutional Sparse Representation and Nonlocal Self-similarity}, author={Weimin Yuana, Yuanyuan Wang, Ruirui Fan, Yuxuan Zhang,
Guangmei Wei, Cai Meng, Xiangzhi Bai}, journal={Computer Vision and Image Understanding}, year={2024}, publisher={Elsevier} }
