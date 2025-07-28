%% Load the data as follow:
%gyroscope reading called: imu_gyr
%accelerometer reading called: imu_acc
%magnetometer reading called: imu_mag
%ground truth quaternion: opt_quat. the quaternion form is: q=[x,y,z,w]
%sampling frequency: sampling_rate


%load the BROAD data trials or any other data her as mintioned earlier:
clear all;
clc;
load('02_undisturbed_slow_rotation_B.mat')
%% processing the data

%delete the non values from the data
nan_rows= any(isnan(opt_quat),2);
opt_quat(nan_rows,:) = [];
imu_acc(nan_rows,:) = [];
imu_gyr(nan_rows,:)=[];
imu_mag(nan_rows,:)=[];
opt_pos(nan_rows,:) = [];

%define the gravitational field value
g=9.7946;

% determine the segment of the data, that will be used to calculate the
% needed variances (the motion during this segment should be static):
st_ind=1;
end_ind=11448;
Std_acc=std(imu_acc(st_ind:end_ind,:));
Std_mag=std(imu_mag(st_ind:end_ind,:));
Std_gyro=std(imu_gyr(st_ind:end_ind,:));
Gyro_cov=diag(Std_gyro).^2;
Mag_cov=diag(Std_mag).^2;
Acc_cov=diag(Std_acc).^2;
% the measurment noise covariance matrix for kalman filter
R=diag([Std_acc,Std_mag]).^2;

a=1
%% Intialization
% define the needed parameters
dt=1/sampling_rate;
seq_len=length(opt_quat(:,1));
% define the start index for the algorithm (usually after the static
% motion)
start_index=11000;
% the predicted quaternion sequence
quat_plus=zeros(seq_len-start_index+1,4);

% intializing the state vector, the covariance matrix and the intial
% quaternion
p0=[1 0 0;0 1 0;0 0 1];
q0=[0 0 0 1];
de0=[0;0;0];
a=1
%% applaying the filter
tic
for i=1:length(quat_plus(:,1))
    j=i+start_index-1;
    [~,quat_plus(i,:),p1]=quaternion_estimation(q0,dt,p0,imu_gyr(j,1),imu_gyr(j,2),imu_gyr(j,3),imu_acc(j,1),imu_acc(j,2),imu_acc(j,3),imu_mag(j,1),imu_mag(j,2),imu_mag(j,3),R,Gyro_cov);
     p0=p1;
     q0=quat_plus(i,:);
end
toc
%% calculate the angles from the quaternion

% the qround truth angles
Ref_Euler=(180/pi).*quat2eul(opt_quat(start_index:seq_len,:));
% the estimated angles 
Euler=(180/pi).*quat2eul([quat_plus(:,4) quat_plus(:,1) quat_plus(:,2) quat_plus(:,3)]);
%calculate the heading_inclination error:
q_error=quatmultiply(opt_quat(start_index:seq_len,:),[quat_plus(:,4) -quat_plus(:,1) -quat_plus(:,2) -quat_plus(:,3)]);
theta=2*acos(sqrt(q_error(:,1).^2+q_error(:,4).^2))*(180/pi);
H_I_error = mean(abs(theta));
a=1
%% plot eulers angles
t=0:dt:dt*(seq_len-start_index);
%plot the heading angle
figure();
plot(t,Euler(:,1),'r');
grid on;
hold on;
plot(t,Ref_Euler(:,1),'b');
legend("estimated","Real");
xlabel("time steps");
title("heading angle");
%plot the heading angle error
figure();
plot(t,(Euler(:,1)-Ref_Euler(:,1)));
grid on;
title("Estimation error of heading angle");
xlabel("time steps");
%calculate the hading error mean and variance
heading_error=abs(Euler(:,1)-Ref_Euler(:,1));
heading_error_mean=mean(heading_error);
heading_error_std=std(heading_error);

%plot the pitch angle
figure();
plot(t,Euler(:,2),'r');
grid on;
hold on;
plot(t,Ref_Euler(:,2),'b');
legend("estimated","Real");
xlabel("time steps");
title("pitch angle");
%plot the pitch angle error
figure();
plot(t,(Euler(:,2)-Ref_Euler(:,2)));
grid on;
title("Estimation error of pitch angle");
xlabel("time steps");
%calculate the pitch error mean and variance
pitch_error=abs(Euler(:,2)-Ref_Euler(:,2));
pitch_error_mean=mean(pitch_error);
pitch_error_std=std(pitch_error);



roll_error=abs(Euler(:,3)-Ref_Euler(:,3));
for l=1:length(roll_error)
    if (roll_error(l)>200)
        roll_error(l)=roll_error(l-1);
    end
end

figure();
plot(t,Euler(:,3),'r');
grid on;
hold on;
plot(t,Ref_Euler(:,3),'b');
legend("estimated","Real");
xlabel("time steps");
title("roll angle");
figure();
plot(t,roll_error);
grid on;
title("Estimation error of roll angle");
xlabel("time steps");
%calculate the roll error mean and variance
roll_error_mean=mean(roll_error);
roll_error_std=std(roll_error);