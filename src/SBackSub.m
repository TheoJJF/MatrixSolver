function x = SBackSub(U,d)
    if iscell(U) == false
        x = BackSub(U,d);
        return
    end

    n = size(U,1);
    x = zeros(n,1);

    x(n) = d(n)/U{n}(n);

    for i = n-1:-1:1
        pivot_recip = 1/U{i}(i);

        for j = i+1:n
            
        end
        

    end
end