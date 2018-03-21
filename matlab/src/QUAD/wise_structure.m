function [ est ] = wise_structure( S_hat, rho, E, options )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Precision Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 4
        input.S = S_hat;
        input.rho = rho;
        input.E = E;
        
        est = wise_structure_main(input, options);
    else
        est = NaN;
        disp('Input not correct!');
    end
end

