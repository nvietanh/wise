function [ output_args ] = wise_structure_print_output( info )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Precision Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Auxiliary function to print iteration output
    iter = length(info.obj);
    display(sprintf('%6i\t%10e\t%10e\t%10e\t%10.2f', iter, info.stepsize(iter), info.norm_grad(iter), info.obj(iter), info.toc(iter)));
end

