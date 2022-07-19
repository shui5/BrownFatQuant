BrownFatQuant Documentation

BrownFatQuant is written in MATLAB to segment and quantify brown adipose tissue (BAT) and white adipose tissue (WAT) using MR images. This work has only been tested on adolescent data acquired using Philips Achieva 3T scanners.


Source code is published on Github https://github.com/shui5/BrownFatQuant

Features:
BrownFatQuant supports segmentation and quantification of BAT and WAT using Gaussian Mixture Model. It fits fat-fraction and T2* signals into the mixture model using expectation–maximization (EM) algorithm by iterating the local maximum likelihood  of them into BAT and WAT. This program requires input of water, fat, fat-fraction and T2* images acquired using the mDixon sequence in Philips or equivalent sequences from other vendors.

Usage:
1) Convert DICOM images to NIFIT (.nii) images.
2) Rename the NIFIT images of water, fat, fat-fraction and T2* to 1-water.nii, 1-fat.nii, 1-ff.nii and 1-t2star.nii.
3) Replace the path of the code and data in run_Brown_Fat_Quant.m to your own path.
4) Run the run_BrownFatQuant.m script (run time varies as the EM algorithm iterating data to fit the mixture model). 
5) Fat-fraction, T2* and volume of BAT and WAT will be reported at the end of the Command window.
6) Masks for BAT and WAT will be generated in the data folder.

Remarks:
The functions folder includes the NIFTI toolbox (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), the GMM toolbox (https://www.mathworks.com/matlabcentral/fileexchange/26184-em-algorithm-for-gaussian-mixture-model-em-gmm) and the denoising toolbox (1. Coupe P, Manjon JV, Gedamu E, Arnold D, Robles M, Collins DL.Robust Rician noise estimation for MR images. Med Image Anal 2010;14:483–493. 2. Manjon JV, Coupe P, Buades A, Louis Collins D, Robles M. New methods for MRI denoising based on sparseness and self-similarity. Med Image Anal 2012;16:18–27.)

For any questions, feedback, suggestions, or critique, please contact Steve Hui <stevehui@jhu.edu>

Should you publish material that made use of BrownFatQuant, please cite the following publication:
Hui SCN, Ko JKL, Zhang T, Shi L, Yeung DKW, Wang D, Chan Q, Chu WCW. Quantification of brown and white adipose tissue based on Gaussian mixture model using water-fat and T2* MRI in adolescents. J Magn Reson Imaging. 2017 Sep;46(3):758-768. doi: 10.1002/jmri.25632. Epub 2017 Jan 16. PMID: 28092409.

Acknowledgements:
This work has been supported by Research Grants Council of the Hong Kong Special Administrative Region; contract grant numbers: 411811; 416712; SEG_CUHK02.
