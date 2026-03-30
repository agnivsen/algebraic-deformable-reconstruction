% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function h = ABFilters1D(filter_name, sigma, hsz);

if ~exist('sigma', 'var') sigma = 1; end;
if ~exist('hsz', 'var') hsz = max(1, 3*sigma); end;

switch filter_name
    
    case 'der_x'
        h = [-1 0 1];
        
    case 'der_y'
        h = [-1 ; 0 ; 1];
        
    case 'gauss_x'
        x = -hsz:hsz;
        h = exp(- (x.^2) / (2*sigma^2)); % / (2*pi*sigma^2);
        h = h / sum(h);
        
    case 'gauss_y'
        x = (-hsz:hsz)';
        h = exp(- (x.^2) / (2*sigma^2)); % / (2*pi*sigma^2);
        h = h / sum(h);

    case 'gauss_der_x'
        x = -hsz:hsz;
        h = x .* exp(- (x.^2) / (2*sigma^2)); % / (2*pi*sigma^2);
        h = h / sum(abs(h));
        
    case 'gauss_der_y'
        x = (-hsz:hsz)';
        h = x .* exp(- (x.^2) / (2*sigma^2)); % / (2*pi*sigma^2);
        h = h / sum(abs(h));        
        
    otherwise
        error('unknown filter');
    
end;
    