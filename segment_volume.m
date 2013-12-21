% Wrapper for the C++ code.
% Johannes Ulén, 2013
function [labelling, lb, e] = segment_volume(unary, include, exclude, conn, conn_weights, options)

%% Parse input
num_regions = length(unary);

% Converts regularization type if zero truncation is used
if (isfield(options,'tau') && isfield(options,'regularization_type'))
	if (options.tau == 0)
		if (~strcmp(options.regularization_type,'scalar'))
			options.regularization_type = 'scalar';
			disp 'Truncation set to zero (tau = 0.0) switching to regularization_type to scalar';
		end
	end
end

% Convert to bool
if isfield(options,'rd')
	options.rd = logical(options.rd);
end

if isfield(options,'improve')
	options.improve = logical(options.improve);
end

if isfield(options,'verbose')
	options.verbose = logical(options.verbose);
end

%% Convert to int32
for iter = 1:length(include)
    include{iter} = int32(include{iter});
		
		% Background must always be included.
		if (include{iter}(end) ~= 1)
			include{iter}(end+1) = 1;
		end
end

for iter = 1:length(exclude)
    exclude{iter} = int32(exclude{iter});
end

%%
% Keeping track on how the unary cost should be assigned
% that is unary_child(5) = 3, means cost \mu_5 - \mu_3
% see the paper for more details.
includes = zeros(num_regions,1);
for inc = 1:length(include)
	constraints = include{inc};
	
	for id = 1:length(constraints)-1;
			child = constraints(id);
			parent = constraints(id+1);
			
			if ( (child < 1) || (parent < 1) ||  ...
					(child > num_regions)  || (parent > num_regions))
				error('ERROR: Trying to include non-existing region.');							
			end
			
			% Not assigned
			if (includes(child) == 0)
				includes(child) = parent;
				
				% See if this is consistent
			else		
				if (includes(child) ~= parent)
					msg1 = sprintf('ERROR: in inclusion constraint %d.', inc);
					msg2 = sprintf(' Region %d is being forced to be contained inside both region %d and %d.', child, parent,  includes(child));
					msg3 = sprintf(' Possible solution: include either region %d inside region %d or vice versa.',  parent,  includes(child));
					error(strcat(msg1,msg2,msg3));
				end
			end
	end
end

% Not contained anywhere set to be contained in background.
includes(includes == 0) = 1;

% Background
includes(1) = 0; 
includes = int32(includes);

%% Check if excluding existing regions
for exc = 1:length(exclude)
    exclude_list = exclude{exc};

    if (any(exclude_list == 1))
       error('ERROR: Cannot exclude the background (region 1)');
    end

    if (any(exclude_list < 0) || any(exclude_list > num_regions))
       error('ERROR: Trying to exclude non-existing region.');	
    end     
end

%% Checks if include and exclude contradict each other.
for exc = 1:length(exclude)
	
	% Excludes in this cell
	exclude_list = exclude{exc};
	
	for inc = 1:length(include)
		
		% Includes in this cell
		include_list = include{inc};
		
		for e_id = 1:length(exclude_list)
			
			% We found this region inside this include list
			id = find(include_list == exclude_list(e_id));
			if ~isempty(id)
				
				all_others = exclude_list;
				all_others(e_id) = [];
				
				for other = all_others
					% all these ids are to be included
					if any(include_list(id:end) == other)
						msg1 = sprintf('ERROR: The inclusion and exclusion lists are contradicting each other.');
						msg2 = sprintf(' Constraints: inclusion %d, exclusion %d. ', inc, exc);
						msg3 = sprintf(' Regions: %d, %d.', include_list(id), other);
						error(strcat(msg1,msg2,msg3));
					end
				end
			end
		end
	end
end


%% Compile (if need be
my_name = mfilename('fullpath');
my_path = fileparts(my_name);

cpp_file = ['cpp' filesep 'segment_volume_mex.cpp'];
out_file = 'segment_volume_mex';
extra_arguments = {['-I"' my_path '"']};
sources = {['cpp' filesep 'maxflow-v3.02.src' filesep 'graph.cpp'], ...
					 ['cpp' filesep 'maxflow-v3.02.src' filesep 'maxflow.cpp'], ...
					 ['cpp' filesep 'QPBO-v1.3.src' filesep 'QPBO.cpp'], ...
					 ['cpp' filesep 'QPBO-v1.3.src' filesep 'QPBO_extra.cpp'], ...
					 ['cpp' filesep 'QPBO-v1.3.src' filesep 'QPBO_maxflow.cpp'], ...
					 ['cpp' filesep 'QPBO-v1.3.src' filesep 'QPBO_postprocessing.cpp']};
				 
compile(cpp_file, out_file, sources, extra_arguments)

%% Solve
[raw_labelling, lb, e] =  ...
	segment_volume_mex(unary, includes, exclude, conn, conn_weights, options);

%% Post process
num_labels = numel(unary{1});
im_size = size(unary{1});

% First region is background: special case.
raw_labelling = reshape(raw_labelling,[im_size num_regions-1]);

if num_regions < 3
	labelling{1} = (raw_labelling == 0);
else
	labelling{1}  = sum(raw_labelling, ndims(raw_labelling)) == 0;
end

for r = 1:num_regions-1
	  inds = num_labels*(r-1)+1:num_labels*r;
    labelling{r+1} = logical( reshape(raw_labelling(inds), im_size));
end

%% Check that the c++ code is working correctly
% Check inclusions
for iter = 1:length(include)
	include_list  = include{iter};
	
	% Removes background
	include_list(include_list == 1) = [];
	
	for i1 = 1:length(include_list)-1
		i2 = i1+1;
		
		A = labelling{include_list(i1)} == 1;
		B = labelling{include_list(i2)} == 1;
		
		if any( (A(:) == 1)  & (B(:) == 0))
			string = sprintf('Inclusion constraint %d is violated  (Region %d not included in %d). \n', ...
												iter, include_list(i2),include_list(i1));
			error('A inclusion constraint is violated');
		end
	end
end

% Check excludes
for iter = 1:length(exclude)
	exclude_list = exclude{iter};
	
	for e1 = 1:length(exclude_list)-1
		for e2 = e1+1:length(exclude_list)
			
			A = labelling{exclude_list(e1)} == 1;
			B = labelling{exclude_list(e2)} == 1;
	
			if any(A(:) & B(:))
				string = sprintf('Exclusion constraint %d is violated  (regions {%d,%d} not excluded). \n', ...
													iter, exclude_list(e1), exclude_list(e2));
				error(string);
			end
		end
	end
end