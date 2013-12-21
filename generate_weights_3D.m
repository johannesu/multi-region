% [I w] = generate_weights_3D(I, delta, D)
%
% I = [dx1 dx2  ...  dxn
%      dy1 dy2  ...  dyn   contains the egdes
%      dz1 dz2  ...  dzn]
%
% I assumed to be symmetric, i.e. only straight lines
%
% delta is the grid spacing (default: 1)
%
% D is the metric (default: eye(3))
%
% NOTE: The functions returns the weights in a different order.
%       Therefore, the sorted I is also returned.
%
% Petter Strandmark 2011
% petter@maths.lth.se
function [I w] = generate_weights_3D(I, delta, D, should_plot, should_print)
	assert(size(I,1)==3);
	if nargin < 2
		delta = 1;
	end
	if nargin < 3
		D = eye(3);
	end
    if nargin < 4
        should_plot = false;
    end
    if nargin < 5
        should_print = should_plot;
    end
	detD = det(D);
    
    %Add voronoi package to path
    path = mfilename('fullpath');
    file = mfilename;
    path = path(1:end-length(file));
    voronoi_path = [path 'sphere_voronoi'];
    addpath(voronoi_path);

	n = size(I,2);
    XYZ = zeros(3, n);
	for i = 1:n
		%Calculate angle
		dx = I(1,i);
		dy = I(2,i);
		dz = I(3,i);
		
		r(i) = sqrt(dx*dx + dy*dy + dz*dz);
		theta(i) = acos(dz/r(i));
		phi(i)   = atan2(dy,dx);
        
        XYZ(1,i) = dx / r(i);
        XYZ(2,i) = dy / r(i);
        XYZ(3,i) = dz / r(i);
	end
	
	%Check symmetry
	for i = 1:n
		dx = I(1,i);
		dy = I(2,i);
		dz = I(3,i);
		corresponding = find( I(1,:)==-dx & I(2,:)==-dy & I(3,:)==-dz);
		if numel(corresponding) ~= 1
			error('Both e and -e need to be present');
		end
	end
	
    %Compute areas of Voronoi diagram
    [ face_num, face ] = sphere_delaunay ( n, XYZ );
    v_xyz = voronoi_vertices ( n, XYZ, face_num, face );
    [ first, list ] = voronoi_polygons ( n, face_num, face );
    list_num = 2 * face_num;
    v_num = face_num;
    A = voronoi_areas ( n, first, list_num, list, XYZ, v_num, v_xyz )';
    
    %Compute weights
    for i=1:n
		dx = I(1,i);
		dy = I(2,i);
		dz = I(3,i);
		e = [dx;dy;dz];
		w(i) = delta^3 * (e'*e)^(3/2) * A(i) * detD / (pi * (e'*D*e)^2);
    end

    if should_plot
        subplot(1,2,1);
        hold off;
        for i = 1:n
            dx = I(1,i);
            dy = I(2,i);
            dz = I(3,i);
            plot3([0 dx], [0,dy], [0 dz], 'k-'); hold on;
            plot3([0 dx], [0,dy], [0 dz], 'k.');
        end
        axis equal
        axis vis3d
        xlim([min(1.5*I(1,:)) max(1.5*I(1,:))]);
        ylim([min(1.5*I(2,:)) max(1.5*I(2,:))]);
        zlim([min(1.5*I(3,:)) max(1.5*I(3,:))]);
%         title(sprintf('N%d',n));
		set(gca,'XTick',[-2 -1 0 1 2]);
		set(gca,'YTick',[-2 -1 0 1 2]);
		set(gca,'ZTick',[-2 -1 0 1 2]);
    
        %Discretize
        npoints=400;
        [THETA PHI] = meshgrid(linspace(0,pi,npoints), linspace(-pi,pi,npoints));
        dtheta = pi/(npoints-1);
        dphi   = 2*pi/(npoints-1);

        bestind = zeros(numel(THETA),1);
        area    = zeros(numel(THETA),1);
        for i = 1:numel(THETA)
            lat1  = pi/2-THETA(i);
            long1 = PHI(i);

            bestdist = inf;
            for j = 1:n
                lat2  = pi/2-theta(j); 
                long2 = phi(j);

                dist = acos( sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(long1-long2) );
                if dist < bestdist
                    bestdist=dist;
                    bestind(i) = j;
                end
            end

            % Area of this element
            area(i) = sin(THETA(i))*dphi*dtheta;
        end

        for i=1:n
            Aapprox(i) = sum( area(bestind==i) );
        end
   
        subplot(1,2,2);
        hold off;
        X=sin(THETA).*cos(PHI);
        Y=sin(THETA).*sin(PHI);
        Z=cos(THETA);
        surf(X,Y,Z,reshape(bestind,size(Z)),'LineStyle','none');
        axis equal
        axis vis3d
        xlim([-1.1 1.1]);
        ylim([-1.1 1.1]);
        zlim([-1.1 1.1]);
%         title(sprintf('N%d',n));
		camva(5.5);
		axis off
        
        A
        Aapprox
    end
    
    if should_print
		for i=1:n
			dx = I(1,i);
			dy = I(2,i);
			dz = I(3,i);
			len2(i) = dx*dx + dy*dy + dz*dz + 0.000001*(100*abs(dx) + 10*abs(dy) + dz);
		end
		[len2 ind] = sort(len2);
		wnew = w(ind);
		A = A(ind); 
		Inew = I(:,ind);
		w = [];
		I = [];
		
        fprintf('\n\n');  
        if n<=10
			fprintf('\\begin{longtable}{rrrll}\n');
			fprintf('\\toprule\n');
			fprintf('$\\Delta x$ & $\\Delta y$ & $\\Delta z$ & Solid angle & $w$\\\\\n');
			fprintf('\\midrule\n');
            for i=1:n
                dx = Inew(1,i);
                dy = Inew(2,i);
                dz = Inew(3,i);
                fprintf('%d & %d & %d & %f & %f \\\\\n',dx,dy,dz,A(i),wnew(i));
            end
		else
			fprintf('\\begin{longtable}{rrrllrrrll}\n');
			fprintf('\\toprule\n');
			fprintf('$\\Delta x$ & $\\Delta y$ & $\\Delta z$ & Solid angle & $w$ &');
			fprintf('$\\Delta x$ & $\\Delta y$ & $\\Delta z$ & Solid angle & $w$\\\\\n');
			fprintf('\\midrule\n');
			for i=1:n/2
                dx = Inew(1,i);
                dy = Inew(2,i);
                dz = Inew(3,i);
                fprintf('%d & %d & %d & %f & %f &',dx,dy,dz,A(i),wnew(i));
				dx = Inew(1,n/2+i);
                dy = Inew(2,n/2+i);
                dz = Inew(3,n/2+i);
				fprintf('%d & %d & %d & %f & %f \\\\\n',dx,dy,dz,A(n/2+i),wnew(n/2+i));
			end
        end
        fprintf('\\bottomrule\n');
        fprintf('\\end{longtable}\n');
        fprintf('\n\n');  
    end
end