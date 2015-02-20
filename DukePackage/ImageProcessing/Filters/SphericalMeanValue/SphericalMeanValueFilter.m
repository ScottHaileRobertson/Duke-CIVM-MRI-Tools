function B_int = SphericalMeanValueFilter(B_delta,mask,radius)
dims = size(B_delta);

[y,x,z] = meshgrid(ceil(-dims(2)/2:dims(2)/2-1),...
    ceil(-dims(1)/2:dims(1)/2-1),...
    ceil(-dims(3)/2:dims(3)/2-1));
clear x_lsp y_lsp z_lsp;
r = sqrt(x.^2 + y.^2 + z.^2);
clear x y z;

% Calculate spherical filter
p = (r<=radius);% & (r>(radius-1));
clear r;
p = p/sum(p(:));
sum_p_1 = sum(p(:))

% Apply spherical mean value filter
smv_filtered = ifftn(fftn(fftshift(p)).*fftn(B_delta));

%Calculate B_delta_prime
B_delta_prime = B_delta - smv_filtered;

% Deconvolve kernel
p(ceil((dims(1)+1)/2),ceil((dims(2)+1)/2),ceil((dims(3)+1)/2))=-1;
p = -1*p;
sum_p_0 = sum(p(:))
B_int = mask.*ifftn(fftn(B_delta_prime)./fftn(fftshift(p)));


% dims = size(B_delta);
% 
% [y,x,z] = meshgrid(ceil(-dims(2)/2:dims(2)/2-1),...
%     ceil(-dims(1)/2:dims(1)/2-1),...
%     ceil(-dims(3)/2:dims(3)/2-1));
% clear x_lsp y_lsp z_lsp;
% r = sqrt(x.^2 + y.^2 + z.^2);
% clear x y z;
% 
% filt = (r<=radius) & (r>(radius-1));
% clear r;
% filt = -1*filt/sum(filt(:));
% filt(ceil((dims(1)+1)/2),ceil((dims(2)+1)/2),ceil((dims(3)+1)/2))=1;
% 
% B_int = ifftn(fftn(filt).*fftn(B_delta));
