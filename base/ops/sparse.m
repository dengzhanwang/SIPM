function C = sparse(A)
    C = unaryOperation(A, @sparse);
end