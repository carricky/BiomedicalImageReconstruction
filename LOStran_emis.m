% remember to run "setup.m" to add the path of PLOT. If you have done so,
% please comment the following two lines.
% PLOTPATH = '../';   % '../' effective only when you are in the "demo" foler 
% run([PLOTPATH,'setup']);  

% image size
imgsiz = [180 180]; 
ntof = 11;
%% loading system matrix
 disp(':: 1. Loading system matrix G ...')
 if ~isvar('GT')
     GTfold = '../GE690TOF2D/GE690_attn_sysmat_180x3.27mm.mat'; % where you store your system matrix G;
     GT = mat_read(GTfold);	% load matrix G
     GTopt.mtype  = 'matlab';
     GTopt.imgsiz = imgsiz;
 end
 if ~isvar('GE')
     GEfold = '../GE690TOF2D/GE690_tof2d_sysmat_11tbin_180x3.27mm.mat';
     GE = mat_read(GEfold);
     GEopt.mtype  = 'tofpet';
     GEopt.imgsiz = imgsiz;
 end

%% simulate PET data
disp(':: 2. Simulating PET data ...')

% load the original image
load('brain_data');
x0 = imresize(x0, imgsiz);
ct = imresize(ct, imgsiz);

% blank scan of transmission
count_blank = 1e5; 
blank = zeros(size(ct));
proj_blank = proj_forw(GT, GTopt, blank);
yi_blank = proj_blank;  % noiseless projection
cs_blank = count_blank / sum(yi_blank(:));
yi_blank = cs_blank * yi_blank;
ni_blank = ones(size(yi_blank))*cs_blank; % normalization factor

% noisy measurements of transmission
count_tran = 1e5; % a total of 100k events
proj_tran  = proj_forw(GT, GTopt, ct); 
% ri = mean(ai(:).*proj(:)) * 0.2 * ones(size(proj));  % 20% uniform background
y0_tran = proj_tran;  % noiseless projection
cs_tran = count_tran / sum(y0_tran(:));
y0_tran = cs_tran * y0_tran;
% ri_tran = cs_tran * ri_tran;
yi_tran = poissrnd(y0_tran); % noisy projection
ni_tran = ones(size(yi_tran))*cs_tran; % normalization factor

% noisy measurements of emission
count_emis = 5e5; % a total of 500k events
proj_emis  = proj_forw(GE, GEopt, x0); 
% nrow = size(proj_emis_init,1)/ntof;
% proj_emis = zeros([nrow, ntof]);
% for i = 1:ntof
%     proj_emis(:,i) = proj_emis_init((i-1)*nrow+1:i*nrow, :);
% end
% ri = mean(ai(:).*proj(:)) * 0.2 * ones(size(proj));  % 20% uniform background
y0_emis = proj_emis; % noiseless projection
cs_emis = count_emis / sum(y0_emis(:));
y0_emis = cs_emis * y0_emis;
% ri_tran = cs_emis * ri_emis;
yi_emis = poissrnd(y0_emis); % noisy projection
ni_emis = ones(size(yi_emis))*cs_emis; % normalization factor


%% reconstruction

disp(':: 3. Image reconstruction ...')

% reconstruction parameters
beta  = 2^-5;
delta = 1e-9;
maxit = 10; 
Gopt.savestep = 10;
xinit = [];
ainit = [];
%the main procedure
disp('reconstruction start!!!')
tic;
x = mlacf(yi_emis, yi_tran, yi_blank, ni_tran, ni_emis, ni_blank, GE, GEopt, GT, GTopt, xinit, ainit, [], [], maxit, beta, ntof, []);
toc;

%% display image
figure, imagesc(reshape(x, imagesiz)); colormap(jet(256));
