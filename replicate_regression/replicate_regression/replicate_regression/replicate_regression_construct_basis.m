function [V, V_reg, W, W_reg] = mr_construct_basis_functions(t,tt,nr,options);

% [V, V_reg, W, W_reg] = mr_construct_basis_functions(t,tt,options);
%
% Build matrices of basis functions


V     = [];  % matrix representing to measured data points (asynchronous)
V_reg = [];  % matrix representing data points to be predicted
             % (time points synchronous and complete)

tmax = max(tt);

switch options.basis
  
  case 'polynomial',
    for it = 1:options.n_comp 
      V     = [V,      (t' /tmax).^it];
      V_reg = [V_reg,  (tt'/tmax).^it];
    end
  
  case 'cos',
    for it = 1:options.n_comp 
      V     = [V,      cos(pi * it * t' /tmax) ];
      V_reg = [V_reg,  cos(pi * it * tt'/tmax) ];
    end
  
  case 'sin',
    for it = 1:options.n_comp 
      V     = [V,      sin(pi * it * t' /tmax) ];
      V_reg = [V_reg,  sin(pi * it * tt'/tmax) ];
    end
  
  case 'sin_half',
    for it = 1:options.n_comp 
      V     = [V,      sin( 1/2 * pi * it * t' /tmax) ];
      V_reg = [V_reg,  sin( 1/2 * pi * it * tt'/tmax) ];
    end
  
  case 'exp',
    for it = 1:options.n_comp 
      V      = [V,      [t'  > 0] .* [1 - exp( -it * t' /tmax) ] ];
      V_reg  = [V_reg,  [tt' > 0] .* [1 - exp( -it * tt'/tmax) ] ];
    end    
  
  case 'sin_horizontal',
    for it = 1:options.n_comp 
      V     = [V,      sin( pi * (it-1/2) * t' /tmax) ];
      V_reg = [V_reg,  sin( pi * (it-1/2) * tt'/tmax) ];
    end
  
  case 'cos+sin',    
    for it = 1:options.n_comp 
      V      = [V,       cos(pi * it * t'/ tmax) ];
      V      = [V,       sin(pi * it * t'/ tmax) ];
      V_reg  = [V_reg,   cos(pi * it * tt'/tmax) ];
      V_reg  = [V_reg,   sin(pi * it * tt'/tmax) ];
    end
    
end

% If necessary, add a basis function representing a smooth jump at the start of the time series
% (last basis function, time scale given by options.t_jump) 

if isfinite(options.t_jump),
    V       = [V,      [ 1 - exp(- t' /options.t_jump)] .* [t' >=0] ];
    V_reg   = [V_reg,  [ 1 - exp(- tt'/options.t_jump)] .* [tt'>=0] ];
end


% Set all basis functions to constant values for negative time values

if options.n_comp * [length(find(tt<0))>0] *  options.constant_before_start,
  iii           = length(find(t<0));
  if find(t>0),
    V(t<0,:)      = repmat(V(iii+1,:),iii,1);
  end
    V(t<0,:)      = 0;
  iii           = length(find(tt<0));
  V_reg(tt<0,:) = repmat(V_reg(iii+1,:),iii,1);
end

% If necessary, add a constant basis function (first basis function)

if options.use_offset,
  V     = [ones(size(t')),V];
  if options.remove_offset,
    V_reg = [zeros(size(tt')), V_reg];
  else
    V_reg = [ones(size(tt')), V_reg];
  end
else,
  V     = [zeros(size(t')),V];
  V_reg = [zeros(size(tt')), V_reg];
end

% ------------------------------------------------------------------------------
% The same, for time derivatives (only if necessary)

W     = [];
W_reg = [];

if options.flag_time_derivative,

switch options.basis
  
  case 'polynomial',
    for it = 1:options.n_comp 
      W     = [W,      it * tmax * (t' /tmax).^[it-1]];
      W_reg = [W_reg,  it * tmax * (tt'/tmax).^[it-1]];
    end
  
  case 'cos',
    for it = 1:options.n_comp 
      W      = [W,     - pi * it/tmax * sin(pi * it * t' /tmax) ];
      W_reg  = [W_reg, - pi * it/tmax * sin(pi * it * tt'/tmax) ];
    end
  
  case 'sin',
    for it = 1:options.n_comp 
      W      = [W,      pi * it/tmax * cos(pi * it * t' /tmax) ];
      W_reg  = [W_reg,  pi * it/tmax * cos(pi * it * tt'/tmax) ];
    end
  
  case 'sin_half',
    for it = 1:options.n_comp 
      W      = [W,      1/2 * pi * it/tmax * cos( 1/2 * pi * it  * t' /tmax)];
      W_reg  = [W_reg,  1/2 * pi * it/tmax * cos( 1/2 * pi * it  * tt'/tmax) ];
    end
  
  case 'exp',
    for it = 1:options.n_comp 
      W      = [W,     [t'  > 0] .* [ it/tmax * exp( - it * t' /tmax)] ];
      W_reg  = [W_reg, [tt' > 0] .* [ it/tmax * exp( - it * tt'/tmax)] ];
    end    
  
  case 'sin_horizontal',
    for it = 1:options.n_comp 
      W      = [W,       pi * (it-1/2)/tmax * cos( pi * (it-1/2) * t' /tmax)];
      W_reg  = [W_reg,   pi * (it-1/2)/tmax * cos( pi * (it-1/2) * tt'/tmax) ];
    end
  
  case 'cos+sin',    
    for it = 1:options.n_comp 
      W      = [W,       - pi * it/tmax * sin(pi * it * t' /tmax)];
      W      = [W,         pi * it/tmax * cos(pi * it * t' /tmax)];
      W_reg  = [W_reg,   - pi * it/tmax * sin(pi * it * tt'/tmax) ];
      W_reg  = [W_reg,     pi * it/tmax * cos(pi * it * tt'/tmax)];
    end
    
end

% If necessary, add a basis function representing a smooth jump at the start of the time series
% (last basis function, time scale given by options.t_jump) 

if isfinite(options.t_jump),
    W       = [W,       1/options.t_jump * exp( - t'  /options.t_jump) ];
    W_reg   = [W_reg,   1/options.t_jump * exp( - tt' /options.t_jump) ];
end


% Set all basis functions to constant values for negative time values

if options.n_comp * [length(find(tt<0))>0] *  options.constant_before_start,
    W(t<0,:)      = 0;
    W_reg(tt<0,:) = 0;
end
  
% If necessary, add a constant basis function (first basis function)

W     = [zeros(size(t')),  W    ];
W_reg = [zeros(size(tt')), W_reg];

end
