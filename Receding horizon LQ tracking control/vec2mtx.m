function mtx = vec2mtx(vec)
    [s1,s2] = size(vec);
    
    mtx = reshape(vec,(s1),[]);
    %mtx = reshape(vec,(s1),[])';
end