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
Std_gyro=std(imu_gyr(st_ind:end_ind,:));
Gyro_cov=diag(Std_gyro).^2;
Acc_cov=diag(Std_acc).^2;
R=Acc_cov;
Q=[Gyro_cov,zeros(3,3);zeros(3,3),Gyro_cov];
%% Intialization
% define the needed parameters
dt=1/sampling_rate;
seq_len=length(opt_quat(:,1));
% define the start index for the algorithm (usually after the static
% motion)
start_index=11000;
pose_plus=zeros(seq_len-start_index+1,6);
Euler=zeros(seq_len-start_index+1,2);

% intializing the state vector and the state vector covariance matrix
p0=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];
pose0=[0 0 9.82  0 0 0];
pose_plus(1,:)=pose0;


%% applying the filter
tic
for i=1:length(pose_plus(:,1))
    j=i+start_index-1;
    [pose_plus(i,:),p1_plus,Euler(i,1),Euler(i,2)]=second_method(pose0,p0,dt,0.9,imu_gyr(j,1),imu_gyr(j,2),imu_gyr(j,3),imu_acc(j,1),imu_acc(j,2),imu_acc(j,3),Q,R);
     p0=p1_plus;
     pose0=pose_plus(i,:);
end
toc

%% calculate the angles
% the predited angles
Euler=(180/pi).*Euler;
% the ground truth angles
Ref_Euler=(180/pi).*quat2eul(opt_quat(start_index:seq_len,:));

%% ploting the results
t=0:dt:dt*(seq_len-start_index);
time=t;
figure();
plot(time,Euler(:,2),'r');
grid on;
hold on;
plot(time,Ref_Euler(:,3),'b');
legend("estimated","Real");
xlabel("time (s)");
ylabel('degree');
title("roll angle");
figure();
plot(time,abs(Euler(:,2)-Ref_Euler(:,3)));
grid on;
title("Estimation error of roll angle");
xlabel("time (s)");
ylabel('degree');
roll_error=abs(Euler(:,2)-Ref_Euler(:,3));
roll_error_mean=mean(roll_error);
roll_error_std=std(roll_error);




figure();
plot(time,Euler(:,1),'r');
grid on;
hold on;
plot(time,Ref_Euler(:,2),'b');
legend("estimated","Real");
xlabel("time(s)");
ylabel('degree');
title("pitch angle");
figure();
plot(t,(Euler(:,1)-Ref_Euler(:,2)));
grid on;
title("Estimation error of pitch angle");
xlabel("time (s)");
ylabel('degree');
pitch_error=abs(Euler(:,1)-Ref_Euler(:,2));
pitch_error_mean=mean(pitch_error);
pitch_error_std=std(pitch_error);
