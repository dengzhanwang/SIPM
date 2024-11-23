% Specify the size of the matrix
n = 5;  % You can change n to any integer >= 3

% Initialize the matrix with zeros
A = zeros(n, n);

% Generate random values for the upper triangle (excluding the last column)
for i = 1:n-1
    for j = i+1:n-1
        A(i, j) = randn();       % Generate random numbers (normal distribution)
        A(j, i) = A(i, j);       % Ensure symmetry
    end
end

% Compute the row sums for rows 1 to n-1
s = sum(A(1:n-1, 1:n-1), 2);

% Adjust the last column to make the row sums zero
A(1:n-1, n) = -s;
A(n, 1:n-1) = A(1:n-1, n)';     % Ensure symmetry for the last row and column

% The diagonal remains zero
% Verify that the diagonal elements are zero
diag_elements = diag(A);

% Verify that the matrix is symmetric
is_symmetric = isequal(A, A');

% Verify that the row sums are zero
row_sums = sum(A, 2);

% Display the results
disp('Generated matrix A:');
disp(A);
disp('Diagonal elements (should be zero):');
disp(diag_elements);
disp('Row sums (should be zero):');
disp(row_sums);
disp(['Is the matrix symmetric? ', num2str(is_symmetric)]);
