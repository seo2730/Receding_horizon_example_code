function mtx = vec2mtx(vec)
    [s1,s2] = size(vec);
    n = sqrt(s2);
    mtx = reshape(vec',n,[]);
    %mtx = reshape(vec,(s1),[])';
end