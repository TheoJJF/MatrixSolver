function [U,d] = SGEPP(A,b)
%{

Perform Sparse Gaussian Elimination with Partial Pivoting (SGEPP).

Parameters: A (Struct/Matrix), b (Vector)
Note: A is assume to be nonsingular
    Error is thrown when the matrix is near/singular.

If A is of the form

A = 

  struct with fields:

    sparse_matrix: {4×1 cell}
                n: 4
                m: 4

A.sparse_matrix = 

  4×1 cell array

    {[dictionary (double --> double) with 2 entries]}
    {[dictionary (double --> double) with 3 entries]}
    {[dictionary (double --> double) with 3 entries]}
    {[dictionary (double --> double) with 2 entries]}

Else, perform standard GEPP

A =

     2    -1     0     0
    -1     2    -1     0
     0    -1     2    -1
     0     0    -1     2


Returns: U (Struct/Matrix), d (Vector)

%}

    if isstruct(A) == false
        [U,d] = GEPP(A,b);
        return
    end

    if A.n ~= size(b,1)
        error("SGEPP:IncompatibleSystem","Incompatible system.");
    end

    if A.n ~= A.m
        warning("SGEPP:NonsquareMatrix","Matrix is not square.");
    end
    
    n = size(A.sparse_matrix,1);

    if A.n ~= n
        error("SGEPP:SingularMatrix","Matrix is singular.");
    end

    
end