function x = Solve(A,b,TOL)
    arguments
        A 
        b
        TOL = 1e-12
    end

    [~,U,~,d] = LU(A,b,TOL);
    x = BackSub(U,d);
end
