function [ta, first_acf_root]=compute_ACF(y,p)
    % ACF - Compute Autocorrelations Through p Lags
    % >> myacf = acf(y,p) 
    %
    % Inputs:
    % y - series to compute acf for, nx1 column vector
    % p - total number of lags, 1x1 integer
    % Output:
    % ta - px1 vector containing autocorrelations
    % first_acf_root - lag at which ACF reaches 0.5
    
    % value series
    [n1, n2] = size(y);
    if (n2~=1)
        error('Input series y must be an nx1 column vector')
    end
    N=n1;
    ybar = mean(y); 
    
    % total number of lags p
    [a1, a2] = size(p);
    if ~((a1==1 & a2==1) & (p<n1))
        error('Input number of lags p must be a 1x1 scalar, and must be less than length of series y')
    end

    % Collect ACFs at each lag i
    ta = zeros(p,1);
    for i = 1:p  
       ta(i) = acf_k(y,i,N,ybar) ; 
    end
    ta=[1; ta];  % mind that ta has size (p+1,1) !

    % 19/11/21 search for first point at which ACF halves -------------- no longer first root of autocorrelation function (the one closest to 0)
	% if r_k[delay-1] r_k[delay]>0 are both positive or negative do nothing: r_k[delay-1]*r_k[delay]>0
	ACF_first_root_found=0;
    first_acf_root=double(-1.0);  % must be double!
    %data_used   % data_used(...,0) is the independent variable
    max_lag=p;
    for lag = 1:max_lag
        if (ACF_first_root_found==0)
            % if ta[lag+1]*ta[lag]>0 are both positive or negative do nothing
            if (abs(ta(lag+1)-0.5)<1.0E-12)  % lag+1 as array indexes start from 1, which is lag 0!
                %disp('Precise zero found!')
                first_acf_root=lag+1.0;  %mind! Only for n_var=1!
                ACF_first_root_found=1;
            else
                if ((ta(lag)-0.5>0.0) && (ta(lag+1)-0.5<0.0))
                    %disp('Approximate zero found!')
                    first_acf_root=double(lag)+0.5; %mind! Only for n_var=1!
                    ACF_first_root_found=1;
                end
            end
        end
    end
    
    % returned first_acf_root is the first point at which ACF reaches 0.5
    % mind that is still expressed as lag, so not a point of the 
    % original independent space

end


% function to compute autocorrelation function values at given lag k
function ta2 = acf_k(y,k,N,ybar)
% ACF_K - Autocorrelation at Lag k
% acf(y,k)
%
% Inputs:
% y - series to compute acf for
% k - which lag to compute acf
% N - no of values in the series (dimension of vector y)
% ybar - mean of values series 

cross_sum = zeros(N-k,1) ;
% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (y(i)-ybar)*(y(i-k)-ybar) ;
end
% Denominator, unscaled variance
yvar = (y-ybar)'*(y-ybar) ;
ta2 = sum(cross_sum) / yvar;

end
