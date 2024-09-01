classdef rubust
   properties
      X
      U
      A
      B
      Q
      R
      W
      gamma
      M
      Lambda
   end
   
   methods
       function obj = rubust(X, U, A, B, Q, R, W, gamma, M, Lambda)
          obj.X = X;
          obj.U = U;
          obj.A = A;
          obj.B = B;
          obj.Q = Q;
          obj.R = R;
          obj.W = W;
          obj.gamma = gamma;
          obj.M = M;
          obj.Lambda = Lambda;
       end
      function cost = cost(obj)
          cost = 0;
          for i = 1:size(obj.U, 1)
              cost = cost + obj.X(i, :)*obj.Q*obj.X(i, :).' + obj.U(i, :)*obj.R*obj.U(i, :).' - obj.gamma^2*obj.W(i, :).*obj.W(i, :);
          end
          cost = obj.X(size(obj.X,1), :)*obj.Q*obj.X(size(obj.X,1), :).';
      end
     

      function obj = updateMandL(obj)
          
          obj.M(size(obj.M, 1), :, :) = obj.Q;

          for i = size(obj.M, 1)-1:-1:1
              
              obj.Lambda(i, :, :) = eye(size(obj.A, 2)) + (obj.B*obj.R^-1*obj.B.' - obj.gamma^-2*eye(size(obj.A, 2)))*squeeze(obj.M(i+1, :, :));
              %test = (obj.B*obj.R^-1*obj.B.' - obj.gamma^-2*eye(size(obj.A, 2)))*squeeze(obj.M(i+1, :, :));
              %disp(test)
              obj.M(i, :, :) = obj.Q + obj.A.'*squeeze(obj.M(i+1, :, :))*squeeze(obj.Lambda(i, :, :))^-1*obj.A;
              
              
          end
      end

      function obj = getU(obj)
          for i = 1:size(obj.U, 1)
              %disp(-obj.R^-1*obj.B.'*squeeze(obj.M(i+1, :, :))*squeeze(obj.Lambda(i, :, :))^-1*obj.A*obj.X(i, :).');
              obj.U(i, :) = -obj.R^-1*obj.B.'*squeeze(obj.M(i+1, :, :))*squeeze(obj.Lambda(i, :, :))^-1*obj.A*obj.X(i, :).';
              
              obj.X(i + 1, :) = obj.A*obj.X(i, :).' + obj.B*obj.U(i, :).' + .5*(rand(size(obj.A, 2), 1)*2-1);
              
              
          end
      end
   end
  
end


