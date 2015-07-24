function phase_unwrapped = unwrap_phase_laplacian(phase_im)
dims_vol=size(phase_im);
[ky,kx,kz] = meshgrid(-dims_vol(2)/2:dims_vol(2)/2-1,...
    -dims_vol(1)/2:dims_vol(1)/2-1,...
    -dims_vol(3)/2:dims_vol(3)/2-1);
kx = fftshift(kx);
ky = fftshift(ky);
kz = fftshift(kz);

k_sq = kx.^2 + ky.^2 + kz.^2;

% Unwrap phase using laplacian.
a = fftn(sin(phase_im));
b = fftn(cos(phase_im));
c = ifftn(k_sq.*a);
d = ifftn(k_sq.*b);

phase_unwrapped = fftn(cos(phase_im).*c - sin(phase_im).*d)./k_sq;
phase_unwrapped(k_sq==0) = 0;
phase_unwrapped = real(ifftn(phase_unwrapped));