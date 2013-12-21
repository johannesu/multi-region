clear all; close all;

%% Simple example
% 6 region (including background) problem.
%
% Submodular constraints:
% ----
% Region 3 forced to be contained inside region 2.
% Region 4 forced to be contained inside region 2.
% Region 5 forced to be contained inside region 3.
% Region 6 forced to be contained inside region 4.
%
% Non-submodular constraints
% -----
% Region 3 and 4 cannot overlap.
% Similar to the heart model in the paper.
% Johannes Ulén 2013. 

%% Settings
% Struct containing all options    

% Output progress for the optimization? Default: false;
options.verbose = true;

% Regularization types
% For region i and	for voxel pair {p,q} we define
% p = unary{i}(p) and q = unary{i}(q).
% There are three different types of regularization implemented
%
% scalar (default)	: w*lambda
% power							: w*lambda * min ( abs(p - q), tau)^beta;
% reciprocal				: w*lambda * ( 1 / ( 1 + beta * min( abs(p - q), tau) );
%
% w is a weighting described in 
% "What metrics can be approximated by geo-cuts, or global optimization of length/area and flux."
% by Vladimir Kolmogorov, and Yuri Boykov. 
% The weight balances up the regularization grids. 
%
% if tau < 0, the truncation is turned off (min ( abs(p - q), tau) -> abs(p-q)).
% 
% If a different regularization is wanted  please modify  cpp/segment_volume_mex.cpp. 
options.regularization_type = 'reciprocal';

% Regularization parameters
options.lambda = 5;
options.beta = 5;
options.tau = -1; % no truncation.

% Solver:
% rd = false. Use lagrangian duality (default)
% rd = true.  Use Roof duality as solver
options.rd = false;

% Maximum number of iterations:
% If rd = false: Maximum supergradient steps.
% If rd = true:  Maximum number of improve steps.
%
% If you set this to 1 and use Lagrangian duality
% you remove, it's equivalent to removing all exclusion constraints.
options.maxiter = 25; 

% If Roof duality is chosen and not all labels are labelled, try 
% to improve the results with the "improve" heuristic?
% This will always result in an equal or better solution but may
% take long time. Number of iterations is controlled by maxiter
options.improve = true;

% Step size rule for lagrangian duality
% "adaptive" (default): is the step sized discussed in the paper.
% it starts at 1 and is cut in half each time the supergradient
% switches sign.
% "oneoverk": the classical 1/k with better theoretical properties
% but much worse convergence in practice.
options.steprule = 'adaptive';
   
%% Define inclusions and exclusions
% Cell array where each cell defines inclusions going from 
% child to parent. e.g 
% include{1} = [4 2 1] 
% Forces region 4 to be contained inside region 2.
% region 2 to be contained inside region 1 
% The order of the regions is crucial.
% N.B recall that region 1 is reserved for background.

% This models the inclusions discussed at the top of this file.
include{1} = [6 4 2 1];
include{2} = [5 3 2 1];

%  Uncomment line below to have zero inclusions
%include = cell(0);

% Each voxel can only be assigned _one_ region 
% of the listed regions in each cell.
% E.g.
% exclude{1} = [1 2 3] 
% Is equivalent to
% exclude{1} = [1 2]
% exclude{2} = [2 3]
% exclude{3} = [1 3]
% 
% The difference is that the latter for the latter
% the code will use 3 times as many lagrangian multipliers.
% The order regions is not important.

% Exclusions discussed at the top of this file.
exclude{1} = [3 4];

% Uncomment line below to have zero exclusions
% exclude = cell(0);

%% Generate a 3D volume
im = generate_data(50);

% The unary cost can be either single (float) or double
im = single(im);

%% Setting up Unary term 
% N.B. The code assumes that unary{1} is background.
% Unary cost is squared distance to a mean
m0 = 1; m1 = 2; m2 = 3; m3 = 3;m4 = 4; m5 = 4;

unary{1} = (im-m0).^2;
unary{2} = (im-m1).^2;
unary{3} = (im-m2).^2;
unary{4} = (im-m3).^2;
unary{5} = (im-m4).^2;
unary{6} = (im-m5).^2;

% Spatial priors (half spaces)
halfway = round(size(im,2)/2);

for pixel = halfway:size(im,2);
    frac =  abs((size(im,2)-pixel)/halfway);
    unary{3}(:,pixel,:) = unary{3}(:,pixel,:)*frac;
end

for pixel = 1:halfway
    frac =  abs(pixel/halfway);
    unary{4}(:,pixel,:) = unary{4}(:,pixel,:)*frac;
end

%% Define connectivity used for regularization
% See generate_weights_3d for more info.
I4 = [1 -1 0 0;
	  0  0 1 -1];
  
I6 = [1 -1 0  0 0  0;
	  0  0 1 -1 0  0;
	  0  0 0  0 1 -1];
    
I10 = [I6(1,:) 1 -1  1 -1;
	   I6(2,:) 1  1 -1 -1;
	   I6(3,:) 0  0  0  0];

% 18 Connectivity with correct weights
I18 = [I10(1,:) I4(1,:)    I4(1,:);
	   I10(2,:) I4(2,:)    I4(2,:);
	   I10(3,:) ones(1,4)  -ones(1,4)];
  
% If the images are anisotropic change this value.    
% see generate_weights_3D.
D = diag(ones(3,1));

% Correct weighting
[conn, conn_weights] = generate_weights_3D(I18, 1, D);


%% Call solver
[labelling, lower_bound, energy] = ...
	segment_volume(unary,include, exclude,  conn, conn_weights, options);
					
%% Show result
figure(1)
clf;
title(sprintf('Energy %g Lower bound: %g \n', energy,lower_bound));
z = round(size(im,3)/2)-1;
plot_regions(labelling, im, z);