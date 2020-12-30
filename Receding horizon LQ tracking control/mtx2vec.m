function vec = mtx2vec(mtx)
    % convert a h by w by d matrix to d by h*w matrix
    [h,w,d] = size(mtx);
    len = h*w*d;
    if isa(mtx, 'gpuArray')
      vec = zeros(len,1, 'gpuArray');
    else
      vec = zeros(len,1);
    end
end
