function [xk_plus,Pk_plus,roll,pitch] = second_method(xk_,Pk_,dt,ca,Gyrox,Gyroy,Gyroz,accx,accy,accz,Q,R)

%xk_ the state vector at time step k-1
%Pk_ the covariance matrix at time step k-1
%dt the sampling period 
%xk the state vector at time step k (prior estiamtion)
%xk_plus the state vector at time step k (posterior estimation)
%Pk the covariance matrix at time step k (prior estiamtion)
%Pk_plus the covariance matrix at time step k (posterior estiamtion)

% calculate the F_(k-1) and G_(k-1) matrices for the state equation:
% xk=F_(k-1)*xk_+G_(k-1)*w
omega=[0 -Gyroz Gyroy;
       Gyroz 0 -Gyrox;
       -Gyroy Gyrox 0];
   
 F_k= cat(1,expm(-dt.*omega),zeros(3,3));
 F_k= cat(2,F_k,cat(1,zeros(3,3),eye(3,3).*ca));
 
 G_k=[0 xk_(3)*dt -xk_(2)*dt 0 0 0;
      -xk_(3)*dt 0 xk_(1)*dt 0 0 0;
      xk_(2)*dt -xk_(1)*dt 0 0 0 0;
      0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1];

% the prior estimate of the state vector and the state covariance matrix
xk=F_k*xk_';
Pk=F_k*Pk_*F_k'+G_k*Q*G_k';
%calculate Kalman gain
H=cat(2,eye(3,3),eye(3,3));
K=Pk*H'*(H*Pk*H'+R)^(-1);
%the posterior estimate of the state vector and the state covariance matrix
y=[accx;accy;accz];
xk_plus=xk+K*(y-H*xk);
Pk_plus=(eye(6,6)-K*H)*Pk*(eye(6,6)-K*H)'+K*R*K';
%estimating the roll and pitch angles

pitch=atan2(xk_plus(2),xk_plus(3));
roll=-asin(xk_plus(1)/9.81);
end

