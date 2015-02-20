% FERMIFILTER
% Applies a fermi filter to an image volume or k-space image if 'fft' flag
% is provided
%
% data = FermiFilter(data, fermi_rolloff, fermi_radius)
%     Applies fermi filter to image volume and returns a filtered image volume
% data = FermiFilter(data, fermi_rolloff, fermi_radius, 'fft')
%     Applies fermi filter to kspace volume and returns a filtered kspace
%     volume.
%
%   Default filter params:
%       -1400 fermi_rolloff=0.9, fermi_radius=0.1
%       -1700 fermi_rolloff=0.15, fermi_radius=0.75
%
function data = FermiFilter(data, fermi_rolloff, fermi_radius, varargin)
% Put data in frequency domain unless fft flag provided
if(nargin<4 | strcmp('fft',varargin{1}))
    data = fftshift(fftn(data));
end

data = data.*FermiFilterGenerator(size(data),fermi_rolloff, fermi_radius);

% Put data back to image domain unless fft flag provided
if(nargin<4 | strcmp('fft',varargin{1}))
    data = ifftn(fftshift(data));
end
% 
% abs_vol = abs(recon_vol);
% abs_vol = abs_vol - min(abs_vol(:));
% abs_vol = abs_vol / max(abs_vol(:));
% fft_vol_abs = fftshift(fftn(abs_vol));
% filt_fft_abs = fft_vol_abs.*fermi_filter;
% filt_vol_abs = ifftn(fftshift(filt_fft_abs));
% figure();imslice(real(filt_vol_abs),'Mag Imag');
% 
% fft_vol_real = fftshift(fftn(real_vol));;
% filt_fft_real = fft_vol_real.*fermi_filter;
% filt_vol_real = ifftn(fftshift(filt_fft_real));
% figure();imslice(real(filt_vol_real),'Filt Real');
% 
% fft_vol_imag = fftshift(fftn(imag_vol));
% filt_fft_imag = fft_vol_imag.*fermi_filter;
% filt_vol_imag = ifftn(fftshift(filt_fft_imag));
% figure();imslice(real(filt_vol_imag),'Filt Imag');