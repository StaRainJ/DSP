function [C0,idx11,idx22] = PreProcedure(X0,Y0,IfPreProcess)
%%  X0::Nx2,  X::3xN

    N0 = size(X0,1);
    C0 = 1:N0;
    idx11 = C0;
%% Data procesing
% Remove repeat data and obvious outliers, But potential inliers would be recalled later!  
if IfPreProcess >= 1
 
    [idxUnique,~] = removeRepeat(X0,Y0); 
    idx11 = idxUnique; 
    N = size(X0,1);
    C0 = idx11;

    if IfPreProcess >= 2
       % Using local consistence verification  
        if N >50
        [C0] = DataInitialization1(X0,Y0,idx11);
        end 
    end

end

idx22 = C0;

