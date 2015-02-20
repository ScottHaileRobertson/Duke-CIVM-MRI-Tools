function [ img ] = BF_2D( img, radius, r_sig, d_sig )
%2D Bilateral Filter
% Ref:
% C. Tomasi and R. Manduchi
% "Bilateral Filtering for Gray and Color Images"
% Proceedings of the 1998 IEEE International Conference on Computer Vision
%
% Inputs:
% img: grayscale image to filter
% radius: filter radius
% r_sig: photometric standard deviation
% d_sig: geometric standard deviation
%
% Output:
% img: filtered version of input "img"

    % Conver to double precision for processing
    img = double(img);

    % Gaussian domain
    [x y] = meshgrid(-radius:radius,-radius:radius);
    D = exp(-(x.^2+y.^2)./(2*d_sig^2));
    D = D(:);
    
    % Mirror the input image to handle edges
    img2 = [img(:, radius+1:-1:2) img img(: ,end-1:-1:end-radius)];
    img2 = [img2(radius+1:-1:2,:); img2; img2(end-1:-1:end-radius,:)];
    
    % Save size before converting to vector
    [x y] = size(img);
    img = img(:);
    
    % Constants to save computation
    q = 2*radius;
    constt = -1/(2*r_sig^2);
    
    for n = 1:length(img)
        
        % Linear to square index conversion
        j = floor((n-1)/x)+1;
        k = n-(j-1)*x;
        
        % Extract central pixel
        I = img2(k:k+q,j:j+q);
        mid = I(radius+1,radius+1);
        I = I(:);
        
        % Compute weights
        F = exp(constt*(I-mid).^2).*D;
        
        % Apply filter
        img(n) = F'*I/sum(F);
        
    end
    
    img = reshape(img,[x y]);

end

