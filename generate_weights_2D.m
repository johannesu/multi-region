% [I w] = generate_weights_2D(I, delta, D)
%
% I = [dx1 dx2  ...  dxn
%      dy1 dy2  ...  dyn]  contains the egdes
%
% I assumed to be symmetric, i.e. only straight lines
%
% delta is the grid spacing (default: 1)
%
% D is the metric (default: eye(2))
%
% NOTE: The functions returns the weights in a different order.
%       Therefore, the sorted I is also returned.
%
% Petter Strandmark 2011
% petter@maths.lth.se
function [I w] = generate_weights_2D(I, delta, D, should_plot, should_print)
	assert(size(I,1)==2);
	if nargin < 2
		delta = 1;
	end
	if nargin < 3
		D = eye(2);
	end
    if nargin < 4
        should_plot = false;
    end
    if nargin < 5
        should_print = false;
    end
	detD = det(D);

	n = size(I,2);
	for i = 1:n
		%Calculate angle
		dx = I(1,i);
		dy = I(2,i);
		
		angles(i) = atan2(dy,dx);
	end
	
	w = nan*zeros(n,1);
	
	[angles ind] = sort(angles);
	I = I(:,ind);
	
	%Check symmetry
	for i = 1:n
		dx = I(1,i);
		dy = I(2,i);
		corresponding = find( I(1,:)==-dx & I(2,:)==-dy);
		if numel(corresponding) ~= 1
			error('Both e and -e need to be present');
		end
		corr(i) = corresponding;
    end
	
    % BK-method 
% 	pos = find(angles >= 0 & angles < pi);
% 	for i = pos
%       %Forward difference
% 		if i == n
% 			da = pi - angles(i);
% 		else
% 			da = angles(i+1) - angles(i);
%         end
% 		dx = I(1,i);
% 		dy = I(2,i);
% 		e = [dx;dy];
% 		w(i) = delta^2 * da * (e'*e) * detD / (2* (e'*D*e)^(3/2));
% 		
% 		%Add for corresponding edge as well
% 		w(corr(i)) = w(i);
%     end
    
    % Partitioning method
    for i=1:n
        if i == n
            next = 2*pi + angles(1);
        else
            next = angles(i+1);
        end
        
        if i==1 
            prev = angles(n) - 2*pi;
        else
            prev = angles(i-1);
        end
        
        %Central difference
        da(i) = (next + angles(i))/2 - (angles(i) + prev)/2;

		dx = I(1,i);
		dy = I(2,i);
		e = [dx;dy];
		w(i) = delta^2 * da(i) * (e'*e) * detD / (2* (e'*D*e)^(3/2));
    end
    
    if should_plot
        clf;
        for i = 1:n
            dx = I(1,i);
            dy = I(2,i);
            plot([0 dx], [0,dy], 'k-'); hold on;
            plot([0 dx], [0,dy], 'k.');
            text(1.2*dx, 1.2*dy, num2str(w(i)));
        end
        xlim([min(1.5*I(1,:)) max(1.5*I(1,:))]);
        ylim([min(1.5*I(2,:)) max(1.5*I(2,:))]);
        axis equal
    end
    
    if should_print
        fprintf('\n\n');
		if n <= 10
			fprintf('\\begin{longtable}{rrc}\n');
			fprintf('\\toprule\n');
			fprintf('$\\Delta x$ & $\\Delta y$ & $w$\\\\\n');
			fprintf('\\midrule\n');
			for i=1:n
				dx = I(1,i);
				dy = I(2,i);
				fprintf('%d & %d  & %f \\\\\n',dx,dy,w(i));
			end
		else
			fprintf('\\begin{tabular}{rrcrrc}\n');
			fprintf('\\toprule\n');
			fprintf('$\\Delta x$ & $\\Delta y$ & $w$ &');
			fprintf('$\\Delta x$ & $\\Delta y$ & $w$\\\\\n');
			fprintf('\\midrule\n');
			for i=1:n/2
				dx = I(1,i);
				dy = I(2,i);
				fprintf('%d & %d  & %f & ',dx,dy,w(i));
				dx = I(1,n/2+i);
				dy = I(2,n/2+i);
				fprintf('%d & %d  & %f \\\\\n',dx,dy,w(n/2+i));
			end
		end
        fprintf('\\bottomrule\n');
        fprintf('\\end{longtable}\n');
        fprintf('\n\n');
    end
end