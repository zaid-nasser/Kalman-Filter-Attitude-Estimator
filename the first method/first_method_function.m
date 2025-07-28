function [qk,qk_plus,Pk] = quaternion_estimation(qk_,dt,Pk_,Gyrox,Gyroy,Gyroz,accx,accy,accz,magx,magy,magz,R,G_cov)

% qk the quaternion at time k
% qk_ the quaternion at time k-1
% dt the time difference between the (k-1) moment and the k moment
%Gyrox,Gyroy,Gyroz the angular velocities from the Gyro sensor

% the prediction stage
qk=expm((dt/2)*Omega_BIG([Gyrox;Gyroy;Gyroz]))*qk_';
%calculate the state matrix at moment k denoted as Fk and the noise matrix
%denoted as G. the state equation: d_ek= Fk*d_e(k-1) + G*W
G=-(dt/2).*eye(3,3);

F=-a2A([Gyrox,Gyroy,Gyroz]);
Fk=eye(size(F))+dt*F;
Q=G*G_cov*G';
% the measurment equation
C=(eye(3,3)+2*a2A(qk(1:3))*a2A(qk(1:3))-2*qk(4)*a2A(qk(1:3)));
% Broad Data value for manetic field
Bn= [-0.336135502772989;15.418399138332108;-40.986348058281713];
g=9.84;
% estimate the z vector
z_estimated=[C*[0;0;g];C*[Bn(1);Bn(2);Bn(3)]];
%the measured z vector
z_measured=[accx accy accz magx magy magz]';
%the error between the predicted z vector and the measured z vector
dz= z_measured-z_estimated;
%the measurment matrix
H=2.*[a2A(C*[0;0;g]); a2A(C*[Bn(1);Bn(2);Bn(3)])];
% calculate the prior covariance matrix for the state vector
Pk_minus=Fk*Pk_*Fk'+Q;
Pk_minus=0.5*(Pk_minus+Pk_minus');
%calculate kalman gain
K=Pk_minus*H'*(H*Pk_minus*H'+R)^(-1);
%calculate the posterior state vector
de_plus=K*dz;
%calculate the posterior quaternion( after correction)
qk_plus=qk+TETA_BIG(qk)*de_plus;
%calculate the posterior covariance matrix for the state vector
Pk=(eye(3,3)-K*H)*Pk_minus*(eye(3,3)-K*H)'+K*R*K';
end

