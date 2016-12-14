% use: [AR, angle] = FiberOrganization('..\Images\Normal_B7928_10_3C_2.png', 1);

function [AR, angle] = FiberOrganization(Image, plotFigure)

    if iscellstr({Image})
        Image = imread(Image);
    end

    % Smooth edges
    I = im2double(Image);
    margin = size(I,2) * 0.1;
    for i=1:margin
        I(i, 1:size(I,2)) = ((i-1) .* I(i, 1:size(I,2))) / margin;
        I(size(I,1)-i+1, 1:size(I,2)) = ((i-1) .* I(size(I,1)-i+1, 1:size(I,2))) / margin;
        I(1:size(I,1), i) = ((i-1) .* I(1:size(I,1), i) / margin);
        I(1:size(I,1), size(I,2)-i+1) = ((i-1) .* I(1:size(I,1), size(I,2)-i+1)) / margin; 
    end

	% Image FFT
	tf=fft2(I(:,:,1));  % If the image is color, only take first channel
	TFF=(fftshift(tf));
	c1=(log(1+abs(TFF)));
    % Normalized grayscale image (between 0 and 1)
	c = mat2gray(c1);   
    
	% Adaptative threshold using Otsu's method
	umb=1.25*graythresh(c);
	bw=im2bw(c,umb);
	
    % Holes filling
    fill=imfill(bw, 'holes');

    % Denoising
    se=strel('disk', 4);
    sinruido=imopen(fill,se);

    % Average filter
    h1=fspecial('average',[16,16]);
    Promediador1=imfilter(sinruido,h1);
    % second filter
    Promediador2=imfilter(Promediador1,h1);
    Promediador_double=mat2gray(Promediador2);
    % Measure properties of image regions
    [L,Ne]=bwlabel(Promediador_double);
    STATS = regionprops(Promediador_double, 'all');
    AR = STATS.MinorAxisLength / STATS.MajorAxisLength;
    angle = STATS.Orientation+90;
    
    if nargin < 2 || plotFigure == 1
        figure;
        subplot(2,3,1),imshow(Image); title('Original image')
        subplot(2,3,2),imshow(c,[]);  title('2D Fourier Transform')
        subplot(2,3,3),imshow(bw,[]); title('Binarization')
        subplot(2,3,4),imshow(bw,[]); title('Holes filling')
        subplot(2,3,5),imshow(sinruido,[]);title('Denoising')
        subplot(2,3,6),imshow(Promediador_double);title('Filtro promediador'); hold on;
        plot([128-100*cosd(STATS.Orientation) 128+100*cosd(STATS.Orientation)] ,[128+100*sind(STATS.Orientation) 128-100*sind(STATS.Orientation)],'r-.');
        xlabel(strcat('Fiber orientation: ', num2str(angle)));
    end
end