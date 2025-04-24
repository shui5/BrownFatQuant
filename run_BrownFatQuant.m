% The following source code is written by Steve Hui and Zhang Teng
% Department of Imaging and Interventional Radiology
% The Chinese University of Hong Kong
%% Start BAT WAT Measurement
close all
%clear
%clc
tic;

% use "\" for Windows, "/" for MacOS
addpath (genpath('/Users/steve/Documents/MATLAB/Code/BrownFatQuant'));
path = '/Users/steve/Documents/MATLAB/Code/BrownFatQuant/sample-data/2015-05-16/Y5667613/mDixon-Quant/S47441/';
%path = '/Users/steve/Documents/MATLAB/Code/BrownFatQuant/sample-data/2015-05-16/Y6049093/mDixon-Quant/';

water   = load_untouch_nii([path '1-water.nii']);
fat     = load_untouch_nii([path '1-fat.nii']);
t2star  = load_untouch_nii([path '1-t2star.nii']);
ff      = load_untouch_nii([path '1-ff.nii']);

pix_dim     = water.hdr.dime.pixdim(2:4); %obtain pixel dimension to calculate volume
fatbinary   = fat; %provide header infomation 

ff.img      = fat.img.*100./(fat.img + water.img + 1e-6); %generate a new fat fraction image
ff.img      = double(ff.img);
t2star.img  = double(t2star.img);

%% Image preprocessing

% % Rician noise estimation and denosing
% verbose     = 0;
% rician      = 1;
% beta        = 1;
% [hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(ff.img, rician, verbose);
% verbose     = 1;
% ff.img      = MRIDenoisingPRINLM(ff.img, hfinal, beta, rician, verbose);
% 
% verbose     = 0;
% rician      = 1;
% beta        = 1;
% [hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(t2star.img, rician, verbose);
% verbose     = 1;
% t2star.img  = MRIDenoisingPRINLM(t2star.img, hfinal, beta, rician, verbose);

%% Create binaries to generate T2* and FF only images
 
out_1       = int16(fat.img);
out_2D      = reshape(out_1,[size(out_1,1) size(out_1,2)*size(out_1,3)]);
threshold   = graythresh(out_2D);
binary_2D   = imbinarize(out_2D, threshold);
out_3D      = reshape(binary_2D,[size(out_1,1) size(out_1,2) size(out_1,3)]);
map         = find(out_3D < 0);
out_3D(map) = 0; %set negative value to zero
fatbinary.img = double(out_3D);

t2star.img  = fatbinary.img.*t2star.img;
ff.img      = fatbinary.img.*ff.img;

save_untouch_nii(fatbinary,[path 'fat_mask.nii']);
save_untouch_nii(t2star,[path 'newt2star.nii']);
save_untouch_nii(ff,[path 'newff.nii']);

%% EM algorithm for Gaussian mixture model (GMM)
cluster = fat; %provide header infomation 
ffimg = reshape(ff.img,[1 numel(ff.img)]);
t2starimg = reshape(t2star.img,[1 numel(t2star.img)]);

% Correct super high T2* values to mean values
map2         = find(t2starimg>44);% scnh
t2starimg(map2) = 0; % scnh
ffimg(map2) = 0; % scnh
nonZeroIndexes1 = t2starimg ~= 0;  % m is your row vector array of numbers.
Meant2star = mean(t2starimg(nonZeroIndexes1));
%t2starimg(map2) =Meant2star;
nonZeroIndexes2 = ffimg ~= 0;  % m is your row vector array of numbers.
Meanff = mean(ffimg(nonZeroIndexes2));
%ffimg(map2) =Meanff;
% end of correction

data = [ffimg;t2starimg];

k = 4; % predefined as 4 (important parameter)
% GMModel = fitgmdist(X,4);
%emresult = emgm(data,k); % emgm function,(data,k), k = num of clusters
em = reshape(emresult,size(ff.img));
cluster.img = em;
save_untouch_nii(cluster, [path 'cluster.nii']);

data( all( ~any( data), 2 ), : ) = []; % removes all rows with all zero
data( :, all( ~any( data ), 1 ) ) = []; % and columns
gkde2(data) % Bivariate Kernel Density Estimation to plot 3D histogram

%% Separate emgm clusters

% To get header infomation 
BATmask = fat; 
BATmaskff = fat; 
BATmaskt2star = fat; 
WATmask = fat; 
WATmaskff = fat; 
WATmaskt2star = fat; 

% find the correct cluster of BAT
clusterimg = reshape(cluster.img, [prod(size(cluster.img)) 1]);
t2starimg = reshape(t2star.img, [prod(size(t2star.img)) 1]);
ffimg = reshape(ff.img, [prod(size(ff.img)) 1]);

t2starMean=[];
ffMean=[];

%iterate clusters
for c = 1:k
 
    ffMean = [ffMean mean(ffimg(clusterimg==c))];
    t2starMean = [t2starMean mean(t2starimg(clusterimg==c))];
    
end

%t2starMean
%ffMean

% Reference values obtained from "Hu HH et al. JMRI 2013 Reson Imaging. 2013 PMID: 23440739"
t2starExpectation=17.9;
ffExpectation= 50;%72.1;
differences=[];
for c=1:k
    differences=[differences,norm([t2starExpectation,ffExpectation]-[t2starMean(c),ffMean(c)],2)];
end

[dMin,I]=min(differences);

cluster.img(find(cluster.img~=I))=0.0;
cluster.img(find(cluster.img==I))=1.0;
BATmask.img = reshape(cluster.img,size(fat.img)); % return to 3D volume
%save_untouch_nii(BATmask, [path 'BATmask-before.nii']);

BATmaskff.img = BATmask.img.*ff.img;
BATmaskt2star.img = BATmask.img.*t2star.img;

save_untouch_nii(BATmask, [path 'BATmask.nii']);
save_untouch_nii(BATmaskff, [path 'BATmask-ff.nii']);
save_untouch_nii(BATmaskt2star, [path 'BATmask-t2star.nii']);

WATmask.img = fatbinary.img - BATmask.img; % Create WAT mask
WATmask.img(find(WATmask.img<=0))=0.0;

WATmaskff.img = WATmask.img.*ff.img;
WATmaskt2star.img = WATmask.img.*t2star.img;

save_untouch_nii(WATmask, [path 'WATmask.nii']);
save_untouch_nii(WATmaskff, [path 'WATmask-ff.nii']);
save_untouch_nii(WATmaskt2star, [path 'WATmask-t2star.nii']);

%% Output FF, T2* and volume of BAT and WAT

% Get mean FF and T2* of BAT
BATffMean = mean(BATmaskff.img(BATmaskff.img~=0));
BATt2starMean = mean(BATmaskt2star.img(BATmaskt2star.img~=0));
disp(['BAT Mean FF = ' num2str(BATffMean) ' %']);
disp(['BAT Mean T2* = ' num2str(BATt2starMean) ' ms']);

% calculate BAT volume
vol_weight_BAT = pix_dim(1)*pix_dim(2)*pix_dim(3)/1000*sum(sum(sum(BATmask.img)));
disp(['BAT volume = ' num2str(vol_weight_BAT) ' mL']);

% Get mean FF and T2* of WAT
WATffMean = mean(WATmaskff.img(WATmaskff.img~=0));
WATt2starMean = mean(WATmaskt2star.img(WATmaskt2star.img~=0));
disp(['WAT Mean FF = ' num2str(WATffMean) ' %']);
disp(['WAT Mean T2* = ' num2str(WATt2starMean) ' ms']);

% calculate WAT volume
vol_weight_WAT = pix_dim(1)*pix_dim(2)*pix_dim(3)/1000*sum(sum(sum(WATmask.img)));
disp(['WAT volume = ' num2str(vol_weight_WAT) ' mL']);

% save results into a table in .txt
Description     = {'BAT ff';'BAT T2*';'BAT vol.';'WAT ff';'WAT T2*';'WAT vol.'};
BATffMean       = round(BATffMean,2);
BATt2starMean   = round(BATt2starMean,2);
vol_weight_BAT  = round(vol_weight_BAT,2);
WATffMean       = round(WATffMean,2);
WATt2starMean   = round(WATt2starMean,2);
vol_weight_WAT  = round(vol_weight_WAT,2);
Mean            = [BATffMean;BATt2starMean;vol_weight_BAT;WATffMean;WATt2starMean;vol_weight_WAT];
Unit            = {'%';'ms';'mL';'%';'ms';'mL'};
T               = table(Description,Mean,Unit);
t1              = datetime('now');
DateString      = datestr( t1 ) ;
filename        = ['tabledata' '_' DateString '.txt'];
writetable(T,filename);
%type tabledata.txt;

toc
