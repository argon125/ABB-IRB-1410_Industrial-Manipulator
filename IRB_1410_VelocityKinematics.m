clc;clear all;close all;
syms th1 th2 th3 th4 th5 th6;
IRB_1410 = [0 90  475 th1; 150 0 0 th2; 600 90 0 th3; 120 -90 720 th4;0 90 0 th5;0 0 85 th6]
symbols = [th1; th2; th3; th4; th5; th6];
zeroTsix= eye(4);
for i=1:6
    zeroTsix = zeroTsix * homoFromDH(IRB_1410(i,1),IRB_1410(i,2),IRB_1410(i,3),IRB_1410(i,4));
end
angles = [30;60;90;30;90;60];
prac_vals = double(subs(zeroTsix,symbols,angles))

posVec = [ zeroTsix(13) zeroTsix(14) zeroTsix(15)]
velJacob = jacobian(posVec,symbols)
velJacob = double(subs(velJacob,symbols,angles))
angVel = [10;20;10;0;-10;0];
endVel = velJacob*angVel

function mat = zRotHomo(ang)
    mat = [cosd(ang) -sind(ang) 0 0; sind(ang) cosd(ang) 0 0;0 0 1 0;0 0 0 1];
end

function mat = xRotHomo(ang)
    mat = [1 0 0 0;0 cosd(ang) -sind(ang) 0 ;0 sind(ang) cosd(ang) 0;0 0 0 1];
end

function homo = homoFromDH(a,alpha,d,theta)
    homo = transl(0,0,d)*zRotHomo(theta)*transl(a,0,0)*xRotHomo(alpha);
end