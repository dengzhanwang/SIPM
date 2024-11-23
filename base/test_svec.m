% Test script for mysvec.m

% Create a symmetric matrix
M = MatCell({
    [1 2 3; 2 4 5; 3 5 6] ,
    [7 8 9; 8 10 11; 9 11 12] });

% Create a cone object
K = Cone({BasicCone('s', 3), BasicCone('s', 3)});


% Call mysvec function
x = mysvec(K, M, 0);

% Print the output
disp(x);

% Call mysmat function

y = mysmat(K, x, 0);

% Print the output
disp(y);

fprintf('Error: %e\n', norm(M - y) / norm(M));