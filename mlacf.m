function [x, out, Phi] = mlacf(yit, pi, Bi, ni_tran, ni_emis, ni_blank, GE, GEopt, GT, GTopt, x0, a0, sit, ri, maxit, beta, ntof, R)
% function [x, out, Phi] = mlacf(yit, pi, Bi, ni_tran, ni_emis, ni_blank, GE, GEopt, GT, GTopt, x0, a0, sijt, ri, maxit, beta, R)
% INPUT:
%   yit     sinogram of TOF-PET emission scan data
%   pi      sinogram of transmission scan data
%   Bi      blank scan data
%   ni      multiplicative factor (normalization, attenuation), can be empty
%   GE      emission scan system matrix
%   GEopt   option set for GE, can be empty if GE is a matlab sparse matrix
%   GT      transmission scan system matrix
%   GTopt   option set for GT, can be empty if GT is a matlab sparse matrix
%   x0      initial image estimate, can be empty
%   a0      initial transmission factor estimate, can be empty
%   sijt    addtive factor (randoms and scatters) of emission scan, can be empty
%   ri      addtive factor (randoms and scatters) of transmission scan, can be empty
%   maxit   maximum iteration number
%   beta    regularization parameter
%   R       regularization operator
%
% OUTPUT
%   x       image estimate in vector
%   out     output
%   Phi     objective function value
%% check inputs
% imgsiz = R.imgsiz;
imgsiz = [180, 180];
numpix = prod(imgsiz);
if isempty(x0)
    x0 = ones(numpix,1);
end

if isempty(a0)
    a0 = ones(size(GT,1),1);
end
[yit, sit, ni_emis] = sino_preprocess(yit, sit, ni_emis);
[pi, ri, ni_tran] = sino_preprocess(pi, ri, ni_tran);

% set Gopt
% TODO sensitity map using the proj_back has problem with the 3-d matrix,
% maybe change it to 2-d matrix.
% GEopt = setGopt(ni_emis, GE, GEopt);
% GTopt = setGopt(ni_tran, GT, GTopt);
if isempty(maxit)
    maxit = 10;
end

% regularization
if isempty(beta) || nargin<8
    beta = 0;
end
if beta>0
%     C1 = sum(abs(R.C),2);
end

% initialization
x    = max(mean(x0(:))*1e-9,x0(:)); 
% x(~Gopt.mask) = 0;
a    = max(mean(a0(:))*1e-9,a0(:)); 
% yeps = mean(yi(:))*1e-9;
% wx   = Gopt.sens;

% output
if nargin>1
    out = []; Phi = [];
end
% out.xest = zeros(length(x(:)), min(maxit,ceil(maxit/Gopt.savestep)+1));
% t1 = tic;

%% iterative loop
for it = 1:maxit   
    printf('the %d iteration', it)
%   yit = reshape(sum(yijt,2),[size(yit,1),size(yit,3)]);
    
% first update the attenuation factor
    proj = proj_forw(GE, GEopt, x);
    qit = ni_emis.*proj;
    qi = reshape(sum(qit,2),[size(qit,1),1]);
    l1 = (yit.*qit)./(repmat(a,[1,ntof]).*qit+sit);
    l2 = repmat((1/ntof)*pi.*Bi./(a.*Bi+ri),[1,ntof]);
    a = a./(Bi+qi).*sum(l1+l2,2); %update 'a', the attenuation factor
    
    
% second update the activity image
    denr = repmat(a,[1,size(qit,2)]).*qit+sit;
    tt = yit./denr;
    tt = tt(:);
    GEf = Gmatrix(GE);
%    right_temp = GEf.*repmat(tt,[1, numpix]);
    right_temp = bsxfun(@times, tt, GE);
    nrow = size(right_temp,1)/ntof;
    right = zeros([nrow, numpix, ntof]);
    for i = 1:ntof
        right(:,:,i) = right_temp((i-1)*nrow+1:i*nrow, :);
    end
    right = sum(right,3);
    right = sum(a.*right,1);
%     right = zeros(numpix, 1);
%     for j = 1:numpix
%         printf('updating the %d pixel',j)
%         right(j,1) = sum(a.*sum(reshape(GE(:,j),[size(GE(:,j),1)/ntof,ntof]).*yit./denr,2)); % right term of the update
%     end
    right = zeros([nrow, numpix, ntof]);
    GE3d = zeros([nrow, numpix, ntof]);
    for i = 1:ntof
        GE3d(:,:,i) = GE((i-1)*nrow+1:i*nrow, :);
    end
    x = x./(sum((repmat(a,[1,ntof]).*GE3d),1))*sum(yit./denr).*right;  
end
    

end


