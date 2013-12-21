% Same as plot_regions but uses contours only
function plot_regions(labels, im, z)
transparency = 0.5;

% Region 1
subplot(1,2,1)
sz = size(im);
emptymask = zeros(sz(1:2));
imagesc(im(:,:,z));
axis equal;
axis off;
colormap gray(256)
hold on

% To align both images and make the gray levels equal
% we overlay a empty mask
showMaskAsOverlay2(transparency,emptymask,'blue',[],0)

subplot(1,2,2)
imagesc(im(:,:,z));
colormap gray(256);
axis equal; axis off;
hold on

bw{1} = labels{2}(:,:,z);
bw{2} = labels{3}(:,:,z) - labels{5}(:,:,z);
bw{3} = labels{4}(:,:,z) - labels{6}(:,:,z);

bw{1}(bw{1} < 0) = 0;
bw{2}(bw{2}  < 0) = 0;
bw{3}(bw{3} < 0) = 0;
bw{4} = labels{5}(:,:,z);
bw{5} = labels{6}(:,:,z);


showMaskAsOverlay2(transparency,bw{1},'blue',[],0)
showMaskAsOverlay2(transparency,bw{2},'red',[],0)
showMaskAsOverlay2(transparency,bw{3},'green',[],0)
showMaskAsOverlay2(transparency,bw{4},[1 0 1],[],0)
showMaskAsOverlay2(transparency,bw{5},[1 1 0],[],0)


% From Matlab file exchange
function varargout = showMaskAsOverlay2(opacity, mask, overlaycolor, varargin)
% Copyright 2010 The MathWorks, Inc.

% Show segmentation (mask) with user-specified transparency/color as overlay on image
%
% Using optional input DELEMASK argument, one can
% easily show multiple segmentation masks on a single image. 
%
% SYNTAX:
%
% SHOWMASKASOVERLAY(OPACITY, MASK, OVERLAYCOLOR)
%     Operates on the image in the current figure, overlays a
%     MASK of opacity OPACITY and of color OVERLAYCOLOR.
%
% SHOWMASKASOVERLAY(OPACITY, MASK, OVERLAYCOLOR, IMG)
%     Takes a handle to an image, or an image itself.
%
% SHOWMASKASOVERLAY(OPACITY, MASK, OVERLAYCOLOR, IMG, DELEMASKS)
%     DELEMASKS is a logical binary, indicating whether existing masks
%     should be deleted before new masks are displayed. Default is TRUE.
%
% SHOWMASKOVERLAY(OPACITY)
%     If an overlayed mask already exists in the current figure,
%     this shorthand command will modify its opacity.
%
% IMGOUT = SHOWMASKASOVERLAY(...)
%     Returns an RGB image of class double, capturing the combined IMG
%     and OVERLAY(s) as image IMGOUT.
%
% [IMGOUT, HNDLS] = SHOWMASKASOVERLAY(...)
%     Also returns a structure of handles to the original image and
%     generated overlays in the current axes.
%
% INPUTS:
%
%     OPACITY       The complement of transparency; a variable on [0,1]
%                   describing how opaque the overlay should be. A
%                   mask of opacity of 0 is 100% transparent. A mask
%                   of opacity 1 is completely solid.
%     MASK          A binary image to be shown on the image of
%                   interest. (Must be the same size as the image operated
%                   on.)
%     OVERLAYCOLOR  A triplet of [R G B] value indicating the color
%                   of the overlay. (Standard "color strings"
%                   like 'r','g','b','m' are supported.) Default
%                   is red.
%     IMG           (Optional) A handle to an image, or an image. By
%                   default, SHOWIMAGEASOVERLAY operates on the image
%                   displayed in the current axes. (If this argument is
%                   omitted, or if the current axes does not contain an
%                   image, an error will be thrown.)
%
%                   Alternatively, IMG may be an image, in which case a new
%                   figure is generated, the image is displayed, and the
%                   overlay is generated on top of it.
%
%     DELEMASKS     Delete previously displayed masks?
%                   This operates at a figure-level. (Default = 1.) 
%
% OUTPUTS:
%
%     HNDLS         A structure containing handles of all images (including
%                   overlays) in the current axes. The structure will have
%                   fields:
%                      Original:   The underlying (non-overlaid) image in
%                                  the parent axes.
%                      Overlays:   All overlays created by
%                                  SHOWMASKASOVERLAY.
%
% EXAMPLES:
% 1)
%                   I = imread('rice.png');
%                   I2 = imtophat(I, ones(15, 15));
%                   I2 = im2bw(I2, graythresh(I2));
%                   I2 = bwareaopen(I2, 5);
%                   figure;
%                   imshow(I);
%                   showMaskAsOverlay(0.4,I2)
%                   title('showMaskAsOverlay')
%
% 2)          
%                   I = imread('rice.png');
%                   AllGrains = imtophat(I, ones(15, 15));
%                   AllGrains = im2bw(AllGrains, graythresh(AllGrains));
%                   AllGrains = bwareaopen(AllGrains, 5);
%                   PeripheralGrains = AllGrains -imclearborder(AllGrains);
%                   InteriorGrains = AllGrains - PeripheralGrains;
%                   figure;
%                   subplot(2,2,1.5)
%                   imshow(I); title('Original')
%                   subplot(2,2,3)
%                   imshow(I)
%                   showMaskAsOverlay(0.4,AllGrains)
%                   title('All grains')
%                   subplot(2,2,4)
%                   imshow(I)
%                   % Note: I set DELEMASKS explicity to 0 here so  
%                   % 'AllGrains' mask is not cleared from figure 
%                   showMaskAsOverlay(0.4,InteriorGrains,[1 1 0],[],0)
%                   showMaskAsOverlay(0.4,PeripheralGrains,'g',[],0)
%                   title('Interior and Peripheral Grains')

% Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% V 1.0 07/05/2007

error(nargchk(1,5,nargin));
if nargin >= 4
    if ~isempty(varargin{1})
        if ishandle(varargin{1})
            imgax = varargin{1};
        else
            figure;
            imshow(varargin{1});
            imgax = imgca;
        end
    else
        imgax = imgca;
    end
    fig = get(imgax,'parent');
    axes(imgax);
else
    fig = gcf;
end

if nargin == 5
    deleMasks = logical(varargin{2});
else
    deleMasks = true;
end


if nargin == 1
    overlay = findall(gcf,'tag','opaqueOverlay');
    if isempty(overlay)
        error('SHOWMASKASOVERLAY: No opaque mask found in current figure.');
    end
    mask = get(overlay,'cdata');
    newmask = max(0,min(1,double(any(mask,3))*opacity));
    set(overlay,'alphadata',newmask);
    figure(fig);
    return
end

% If the user doesn't specify the color, use red.
DEFAULT_COLOR = [1 0 0];
if nargin < 3
    overlaycolor = DEFAULT_COLOR;
elseif ischar(overlaycolor)
    switch overlaycolor
        case {'y','yellow'}
            overlaycolor = [1 1 0];
        case {'m','magenta'}
            overlaycolor = [1 0 1];
        case {'c','cyan'}
            overlaycolor = [0 1 1];
        case {'r','red'}
            overlaycolor = [1 0 0];
        case {'g','green'}
            overlaycolor = [0 1 0];
        case {'b','blue'}
            overlaycolor = [0 0 1];
        case {'w','white'}
            overlaycolor = [1 1 1];
        case {'k','black'}
            overlaycolor = [0 0 0];
        otherwise
            disp('Unrecognized color specifier; using default.');
            overlaycolor = DEFAULT_COLOR;
    end
end

figure(fig);
tmp = imhandles(fig);
if isempty(tmp)
    error('There doesn''t appear to be an image in the current figure.');
end
try
    a = imattributes(tmp(1));
catch %#ok
    error('There doesn''t appear to be an image in the current figure.');
end
imsz = [str2num(a{2,2}),str2num(a{1,2})]; %#ok

if ~isequal(imsz,size(mask(:,:,1)))
    error('Size mismatch');
end
if deleMasks
    delete(findall(fig,'tag','opaqueOverlay'))
end

overlaycolor = im2double(overlaycolor);
% Ensure that mask is logical
mask = logical(mask);

if size(mask,3) == 1
    newmaskR = zeros(imsz);
    newmaskG = newmaskR;
    newmaskB = newmaskR;
    %Note: I timed this with logical indexing (as currently implemented),
    %with FIND, and with logical indexing after converting the mask to type
    %logical. All three were roughly equivalent in terms of performance.
    newmaskR(mask) = overlaycolor(1);
    newmaskG(mask) = overlaycolor(2);
    newmaskB(mask) = overlaycolor(3);
elseif size(mask,3) == 3
    newmaskR = mask(:,:,1);
    newmaskG = mask(:,:,2);
    newmaskB = mask(:,:,3);
else
    beep;
    disp('Unsupported masktype in showImageAsOverlay.');
    return
end

newmask = cat(3,newmaskR,newmaskG,newmaskB);

hold on;
h = imshow(newmask);
try
    set(h,'alphadata',double(mask)*opacity,'tag','opaqueOverlay');
catch %#ok
    set(h,'alphadata',opacity,'tag','opaqueOverlay');
end
if nargout > 0
    varargout{1} = imhandles(imgca);
end
if nargout > 1
    varargout{2} = getframe;
    varargout{2} = varargout{2}.cdata;
end