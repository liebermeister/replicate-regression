function [t,x,sigma,r,x_clean,t_true,x_true] = demo_replicate_regression_create_data_multiple(method)


% ---------------------------------
% Create artificial data

randn('state',1);

noise_level = 0.1;

%% 16 copies of the same gene
switch method,
  case 'identical',
    t_true  = 0:1:30;
    x_true  = repmat(2 * [0.1*t_true./(1+0.1*t_true)] .* exp(-0.1*t_true/2), 16,1);
    offset2 =  0.1*[x_true];
    offset1 = -0.5*[x_true];
    offset3 =  0.5*[x_true];
  case 'different',
    t_true  = 0:5:150;
    patterns = diag(1./[1:5]) * sin(1/5*pi/4*[0.1:0.1:0.5]' * t_true);
    x_true   = 1 + 0.4*randn(16,5) * patterns;
    offset1  = 0.2*randn(16,5) * patterns;
    offset2  = 0.2*randn(16,5) * patterns;
    offset3  = 0.2*randn(16,5) * patterns;
  otherwise error('');
end

ind1 = 5:6:20;
ind2 = 10:6:25;
ind3 = 18:6:30;

t1 = t_true(ind1); 
t2 = t_true(ind2);
t3 = t_true(ind3);

x1_true_all = x_true + offset1;
x2_true_all = x_true + offset2;
x3_true_all = x_true + offset3;

x1_true = x_true(:,ind1) + offset1(:,ind1);
x2_true = x_true(:,ind2) + offset2(:,ind2);
x3_true = x_true(:,ind3) + offset3(:,ind3);

x1 = x1_true + noise_level * randn(size(x1_true));
x2 = x2_true + noise_level * randn(size(x2_true));
x3 = x3_true + noise_level * randn(size(x3_true));

x1(x1< 0.1) = 0.1;
x2(x2< 0.1) = 0.1;
x3(x3< 0.1) = 0.1;

sigma1 = noise_level * ones(size(x1));
sigma2 = noise_level * ones(size(x2));
sigma3 = noise_level * ones(size(x3));

t     = [t1, t2, t3];
x     = [x1, x2, x3];
x_clean = [x1_true, x2_true, x3_true];
sigma = [sigma1, sigma2, sigma3];
r     = [ones(size(t1)), 2*ones(size(t2)), 3*ones(size(t3))];

[t,order] = sort(t);
x         = x(:,order);
x_clean   = x_clean(:,order);
sigma     = sigma(:,order);
r         = r(order);
