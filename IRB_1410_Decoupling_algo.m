clc;clear all;close all;

syms t1 t2 t3 t4 t5 t6

d=0;
T1=dh_param(30,90,0,0.475);
T2=dh_param(60,0,0.15,0);
T3=dh_param(30,90,0.6,0);
T4=dh_param(60,-90,0.12,0.72);
T5=dh_param(30,90,0,0);
T6=dh_param(60,0,0,0.085+d); %d is the gripper length


final_matrix = T1*T2*T3*T4*T5*T6;
T03=T1*T2*T3;
T36=T4*T5*T6;
rot_mat=final_matrix(1:3,1:3);
pos_mat=final_matrix(1:3,end); 
a=pos_mat(1,1);
b=pos_mat(2,1);
c=pos_mat(3,1);

theta1=atan2d(b,a)+10 

D=(a^2 +b^2 +(c-0.475)^2 - 0.6^2 -0.12^2)/(2*0.6*0.12);
theta3=atan2d(abs(sqrt(1-D^2)),D) +10
theta2=atan2d((c-0.475)^2,(a^2 +b^2))-atan2d(0.12*sind(theta3),(0.6+0.12*cosd(theta3)))+10
inv_0R3=inv((T1*T2*T3));
RHS=inv_0R3(1:3,1:3)*rot_mat;

theta4=atan2d(RHS(2,3),RHS(1,3))
theta5=acosd(RHS(3,3))
%theta6=-1*atan2d(RHS(3,2),RHS(3,1))
theta6=asind(RHS(3,2)/sind(theta5))  
%heta6=-1*acosd(RHS(3,1)/sind(theta5)) 

function dh=dh_param(theta,alpha,a,d)
    dh=transl(0,0,d)*trotz(theta,'deg')*transl(a,0,0)*trotx(alpha,'deg'); 
end