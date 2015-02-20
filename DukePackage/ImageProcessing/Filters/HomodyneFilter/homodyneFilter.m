function filt_im = homodyneFilter(im, sigma)
% Create filter in frequency domain
dims = size(im);
lsp_ky = linspace(-1,1,dims(1)+1); lsp_ky = lsp_ky(1:end-1);
lsp_kx = linspace(-1,1,dims(2)+1); lsp_kx = lsp_kx(1:end-1);
lsp_kz = linspace(-1,1,dims(3)+1); lsp_kz = lsp_kz(1:end-1);
[kx,ky,kz] = meshgrid(lsp_kx, lsp_ky, lsp_kz);
gauss_filt = 1/((sigma*sqrt(2*pi))^3)*exp(-(kx.^2 + ky.^2 + kz.^2)/(2*sigma^2));
gauss_filt = gauss_filt/max(gauss_filt(:));

% Apply low pass filter
low_pass_im = ifftn(fftshift(gauss_filt.*fftshift(fftn(im))));

% Calculate difference
filt_im = abs(im) - abs(low_pass_im);
