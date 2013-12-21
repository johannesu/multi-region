
%% Generating synthetic data
function im= generate_data(sz)
randn('state',0);

% Standard deviation
sigma = 0.1;

% Means
m0 = 1; 
m1 = 2; 
m2 = 3; 
m3 = 3; 
m4 = 4; 
m5 = 4;

im = randn(sz,sz,sz)*sigma+m0;

r = sz/3;
[x, y, z] = ndgrid(1:sz,1:sz,1:sz);
	
log = (x(:) - sz/2).^2 + (y(:) - sz/2).^2 + (z(:) - sz/2).^2 <= r^2;
im(log) = randn(sum(log),1)*sigma+m1;

r = r/2.2;
log = (x(:) - sz/2).^2 + (y(:) - sz/2 + r).^2 + (z(:) - sz/2).^2 <= r^2;
im(log) = randn(sum(log),1)*sigma+m2;

log = (x(:) - sz/2).^2 + (y(:) - sz/2 - r).^2 + (z(:) - sz/2).^2 <= r^2;
im(log) = randn(sum(log),1)*sigma*2+m3;

r = r/2.2;
log = (x(:) - sz/2).^2 + (y(:) - sz/2 + 2.2*r).^2 + (z(:) - sz/2).^2 <= r^2;
im(log) = randn(sum(log),1)*sigma+m4;

log = (x(:) - sz/2).^2 + (y(:) - sz/2 - 2.2*r).^2 + (z(:) - sz/2).^2 <= r^2;
im(log) = randn(sum(log),1)*sigma+m5;

im = im(:,:,round(sz/2)-5:round(sz/2)+5);