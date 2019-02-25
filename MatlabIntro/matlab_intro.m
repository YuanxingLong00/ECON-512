% Econ 627 Lab 1. Introduction to Matlab. 
%#ok<*ASGLU>
%#ok<*CLSCR>

%% Part1: Basics

clear all 

%--------------------------------------------------------------------------
% Matrix creation
%--------------------------------------------------------------------------

V = [1, 3, 4, 2];           % 1x4 vector. Note that ';' suppresses output
                            % e.g. a = 1, b = 2; c = 3, d = 0;
A = [  3,  1, 42; ...       % 3x3 matrix. Note that ellipsis (...) continues 
     100, -4,  3; ...       % a statement to the next line
       0,  7,  1];
   
B = [0.1, 0.1, 2; ...  
     0.5,  -1, 3; ...
       0,  25, 0];
   
T  = 1 : 10;                % 1, 2, 3, ..., 9, 10
TH = 1 : .5 : 10;           % 1, 1.5, 2, ..., 9.5, 10

%--------------------------------------------------------------------------
% Special matrices
%--------------------------------------------------------------------------

Z = zeros(3, 4);
O = ones(5, 2);
E = eye(3);
N = nan(4, 3);
T = true(1, 4); F = false(4, 1);
I = inf(5);
R = rand(4, 3);             % Random matrix: [R]_{ij}~U(0,1)

%--------------------------------------------------------------------------
% Access elements of a matrix/vector
%--------------------------------------------------------------------------

x1 = V(2);                  % Selects [V]_2
x2 = A(2, 2);               % Selects [A]_{2,2}
x3 = A(8);                  % Selects [A]_{2,3}; matrices are stored in comlumn-major order

%--------------------------------------------------------------------------
% Select specific elements of a matrix/vector
%--------------------------------------------------------------------------

x3 = A(1, :);               % Selects row #1
x4 = A(:, 2);               % Selects column #2
x5 = A([1, 3], 1:2);        % Selects rows #1, #3 and columns from #1 to #2 
x6 = TH(~(TH >= 5));        % Selects numbers less than 5. Same as TH(TH < 5)
x7 = TH(round(TH) == TH);	% Selects integers

%--------------------------------------------------------------------------
% Matrix/vector operations
%--------------------------------------------------------------------------

[n, m] = size(A);           % Dims of a matrix
[~, m] = size(A);           % # of columns; put '~' to ignore an argument
n = size(A);                % # of rows and columns
n = numel(A);               % # of elements

VT  = V'; AT = A';          % Transpose
y1  = A * 5;                % Scalar multiplication
y2  = A / 5;                % Scalar division
y3  = A + E;                % Addition
y4  = A - E;                % Subtraction
y5  = A * B;                % Matrix multiplication
y6  = A \ B;                % Matrix left division:  A \ B = inv(A) * B
y7  = A / B;                % Matrix right division: A / B = A * inv(B)
y8  = A .* B;               % Element-wise product (AKA Hadamard/Schur product)
y9  = B ./ A;               % Element-wise division
y10 = A .^ B;               % Element-wise exponentiation
y11 = 2 .^ A;               % ...???
y12 = kron(R, A);           % Kronecker product

C1 = [A; B]; C2 = [A B];    % Concat. matrices
y13 = repmat(V, 2, 5);      % Creates a big matrix consisting of copies of V
y14 = inv(A);               % Matrix inverse 
y15 = trace(A);             % Trace 
y16 = det(A);               % Determinant
y17 = rank(B);              % Rank of a matrix
y18 = diag(A);              % Selects diagonal elements
y19 = diag(V);              % Create a diagonal matrix with y18 on the main diagonal

[v, e] = eig(A);            % Eigendecomposition of a matrix
[s, v, d] = svd(A);         % Singular value decomposition
c = chol(R' * R);           % Cholesky decomposition

%--------------------------------------------------------------------------
% Stats
%--------------------------------------------------------------------------

VS  = sum(V);               % Sum of all elements
VM  = mean(V);              % Sample mean
VV  = var(V);               % Sample variance
VSD = std(V);               % Sample std
[VM, I1] = min(V);          % Minimum elements with its index

RS  = sum(R, 1);            % Column-wise sum of all elements
RM  = mean(R, 2);           % Rown-wise sample mean
RV  = var(R, 0, 1);         % Column-wise sample variance
RSD = std(R, 0, 2);         % Row-wise sample std
[RM, I2] = min(R, [], 2);   % Column-wise minimums of a matrix with their indices
RC1 = cov(R);               % Covariance matrix (column-wise)
RC2 = corr(R);              % Correlation matrix (column-wise)

%--------------------------------------------------------------------------
% Other operations
%--------------------------------------------------------------------------

z1 = any(V < 0); z2 = any(A < 0, 1);
z3 = all(V > 0); z4 = all(A > 0, 2);
z5 = isempty([]);
z6 = find(R(:, 1) > 0.5); 
z7 = isnan([nan, 1, 2]);
z8 = isinf([-inf, 4, 3]);

%--------------------------------------------------------------------------
% Cell arrays
%--------------------------------------------------------------------------

% Cell arrays creation
C1 = cell(3, 4);
C2 = {A, B; R, V};

% Access elements of a cell array
w1 = C2(2, 1);              % Returns {{C2}_{2,1}} 
w2 = C2{2, 1};              % Returns {C2}_{2,1} - actual content of a cell

%--------------------------------------------------------------------------
% Structs
%--------------------------------------------------------------------------

struct.A = 'A';
struct.B = 1;
struct.C = [1; -1; 0];
struct.B = struct.B + 1;

%--------------------------------------------------------------------------
% Anonymous functions
%--------------------------------------------------------------------------

f = @(x, y) exp(x + y);
disp(f(1, 2));

% Useful examples
w3 = bsxfun(@times, R, R(:, 1));    % Element-wise multiplication along the first dim.
                                    % Same as R .* repmat(R(:, 1), 1, 3)
w4 = bsxfun(@rdivide, R, R(1, :));  % Element-wise division along the second dim.
                                    % Same as R ./ repmat(R(1, :), 4, 1)
SA(1).a = rand(3, 6);
SA(2).a = magic(12);
SA(3).a = ones(5, 10);
w5 = arrayfun(@(x) numel(x.a), SA);
                                     
%% Part 2: Control Flow

clear all

%--------------------------------------------------------------------------
% Conditional Control - if, else
%--------------------------------------------------------------------------

a = randi(100, 1);

if a < 30
    disp('small')
elseif a < 80
    disp('medium')
else
    disp('large')
end

%--------------------------------------------------------------------------
% Conditional Control - switch
%--------------------------------------------------------------------------

[dayNum, dayString] = weekday(date, 'long', 'en_US');

switch dayString
   case 'Monday'
      disp('Start of the work week')
   case 'Tuesday'
      disp('Day 2')
   case 'Wednesday'
      disp('Day 3')
   case 'Thursday'
      disp('Day 4')
   case 'Friday'
      disp('Last day of the work week')
   otherwise
      disp('Weekend!')
end

%--------------------------------------------------------------------------
% Loop Control - for
%--------------------------------------------------------------------------

X = zeros(100, 10);

for i = 1:2:100
    for j = 1:10
        X(i, j) = 1 / (i + j);
    end
end

%--------------------------------------------------------------------------
% Loop Control - while
% Using interval bisection to find a zero of a polynomial
%--------------------------------------------------------------------------

a = 0; fa = -inf;
b = 3; fb = inf;
f = @(x) (x ^ 3 - 2 * x - 5);

while b - a > eps * b
   x = (a + b) / 2; fx = f(x);
   if sign(fx) == sign(fa)
      a = x; fa = fx;
   else
      b = x; fb = fx;
   end
end

disp(x)

%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------

doc plot
doc line
doc bar
doc scatter
doc hist

%% Part 3 : Some Econometric Applications
% Reference: Bruce Hansen's textbook -
% http://www.ssc.wisc.edu/~bhansen/econometrics/Econometrics.pdf

%--------------------------------------------------------------------------
% Application 1. Least squares estimator and standard errors.
% Reference: (Hansen Chapter 4, p. 88)
%--------------------------------------------------------------------------

clear all

% Reseed the global random number generator
rng(135)

beta = [1; 0.5; -0.7; 1]; % parameters.
n = 10000;                 % # of obs.
k = numel(beta);           % # of parameters.

% Generate random data
X1 = exprnd(2, n, 1);
X2 = mvnrnd([0; 0], [1 0.2; 0.2 1], n);

X = [X1 X2]; C = ones(n, 1);
U = normrnd(0, 2, n, 1);
Y = [C X] * beta + (3 * U);

params.alpha = 0.05;
params.var_type = '';

[b_hat1, se1, stats1] = ols(Y, X, params);

% Simulate heteroskedasticity
Y = [C X] * beta + (1 + X1) .* U;

% Homoskedastic standard errors
params.var_type = '';
[b_hat2, se2, stats2] = ols(Y, X, params);

% White standard errors
params.var_type = 'white';
[b_hat3, se3, stats3] = ols(Y, X, params);

%--------------------------------------------------------------------------
% Application 2. Nonlinear least squares.
% Reference: (Hansen, Chapter 12, p. 279)
%--------------------------------------------------------------------------

clear all

% Reseed the global random number generator
rng(135)

beta = [1; 2; -1]; % parameters
n = 10000;                 % # of obs.
k = numel(beta);  
S = [ 1.0  0.2 -0.1; ...
      0.2  1.0  0.1; ...
     -0.1  0.1  1.0];

X = mvnrnd(zeros(3, 1), S, n);
U = normrnd(0, 1, n, 1);

% Logistic link regression
m_fn = @(X, b) (1 + exp(-X * b)) .^ -1;
Y = m_fn(X, beta) + U;

% Scatter plot
scatter(X * beta, Y);
title('Scatter plot')
xlabel('$X_i^\top\beta$', 'Interpreter', 'latex')
ylabel('$Y_i$', 'Interpreter', 'latex')

% Estimation
params.alpha = 0.05;
params.solver = 'fminunc';
params.solver_options = ...
    optimset('Display', 'iter', 'MaxFunEvals', 1000, ...
             'TolFun', 10e-8, 'TolX', 10e-8);

beta0 = zeros(k, 1);
[b_hat1, se1, stats1] = nls(m_fn, Y, X, beta0, params);

% Try another solver
params.solver = 'fminsearch';

[b_hat2, se2, stats2] = nls(m_fn, Y, X, beta0, params);
