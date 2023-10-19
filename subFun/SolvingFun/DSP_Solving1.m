function [X, iter, min_cost,E] = DSP_Solving1(grad, proj,U_E, Xinit,Einit, L, opts, calc_F)  


    min_cost = inf;
    ttemp = [];
     Linv = 1/L;    
    lambdaLiv = opts.lambda*Linv; % lambda/L;
    
%% MAIN LOOP   
    opts_proj = opts;
    opts_proj.lambda = lambdaLiv; % for solving x
    
    opts_proj2 = opts;
    opts_proj2.lambda = opts.mu; % for solving e
    
    x_old = Xinit;
    y_old = Xinit;
    E_old = Einit; 
    t_old = 1;
    iter = 0;
   %% ====================================================================
   % proj = L1 solution------solving  X = argmin_X 0.5*||X - U|| + lambda||X||_1 
   % ------>                          X = feval(proj, U, lambda)  
   % grad = M'*M*X + M'*E;
   % U_E =  -M*X;
   %% ====================================================================
    while  iter < opts.max_iter
        
        iter = iter + 1;
        PL = y_old - Linv*feval(grad, y_old,E_old);
        x_new = feval(proj, PL, opts_proj);
        x_new = SphereProj(x_new);
        DX    = feval(U_E, x_new);               

 %%
        E_new = feval(proj, DX, opts_proj2); 

        t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
        
        dx = (x_new - x_old);
        
        y_new = x_new + (t_old - 1)/t_new * dx;

         y_new = SphereProj(y_new);
         
%% check stop criteria
            e = norm1(dx)/numel(x_new);
%             e = norm2_cols(x_new - x_old);
       
        
          if iter> 100 && mod(iter,20) == 0 && opts_proj2.lambda >=0.02
            opts_proj2.lambda = opts_proj2.lambda.*0.98;
            if e > 1e-4 
                opts_proj2.lambda = opts_proj2.lambda.*0.98;
            end            
         end


%        if ifplot ==1 && mod(iter,10) == 0 %|| iter==opts.max_iter ||e < opts.tol
%         aaaa0 = opts_proj2.lambda;
%        aaaa=  [DX,E_new,DX-E_new];
%         figure(33), plot(DX,'b.'), ylabel('MX'), hold on;
%          line([0,length(DX)],[-aaaa0,-aaaa0],'Color',[.8 .1 .1]); 
%         line([0,length(DX)],[aaaa0,aaaa0],'Color',[.8 .1 .1]); hold off
%         set(gcf,'position',[100 400 500 500]);  
%         drawnow;         
%        end
  
        if e < opts.tol
            break;
        end
        
%         CC = feval(calc_F, x_old, E_old,opts.lambda,opts_proj2.lambda); % 
  %% update
%        ttemp = [ttemp,x_old];
        x_old = x_new; 
        t_old = t_new;
        y_old = y_new;
        E_old = E_new;
    

    end
    
      X = x_old;
      E = E_old; 
    min_cost = feval(calc_F, X, E, opts.lambda,opts_proj2.lambda);
    
 %% plot results  
%  if ifplot==1
%     aaaa0 = opts_proj2.lambda;
%        aaaa=  [DX,E_new,DX-E_new];
%         figure(33), plot(DX,'b.'), ylabel('MX'), hold on;
%          line([0,length(DX)],[-aaaa0,-aaaa0],'Color',[.8 .1 .1]); 
%         line([0,length(DX)],[aaaa0,aaaa0],'Color',[.8 .1 .1]); hold off
%         set(gcf,'position',[100 300 300 300]);  
%         drawnow; 
        
%         figure (44),plot(DX-E_new,'b.'), ylabel('MX-E');
%         set(gcf,'position',[100 100 300 300]);  
%         drawnow; 
%  end
 
end
