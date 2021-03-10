clc;clear all;close all;
syms th1 th2 th3 th4 th5 th6;
IRB_1410 = [0 90  475 th1; 150 0 0 th2; 600 90 0 th3; 120 -90 720 th4;0 90 0 th5;0 0 85 th6]
symbols = [th1; th2; th3; th4; th5; th6];
zeroTsix= eye(4);
for i=1:6
    zeroTsix = zeroTsix * homoFromDH(IRB_1410(i,1),IRB_1410(i,2),IRB_1410(i,3),IRB_1410(i,4));
end
zeroTsix
angles = [30;60;90;30;90;60];
prac_vals = double(subs(zeroTsix,symbols,angles))
func = reshape(zeroTsix-prac_vals,[16,1])

assumption = [30 ;90; 60; 90; 30; 45];
%assumption = [0 ;90; 0; -90; 60; 30];
tol = 0.001;numIter = 100;
jac = jacobian(reshape(zeroTsix,[16,1]),symbols);
for i=1:numIter
    hist = assumption;
    mat = double(subs(jac, symbols,assumption));
    func_val = double(subs(func,symbols,assumption));
    assumption = assumption - pinv(mat)*func_val;
    if(norm(hist-assumption)<tol)
        break
    end
end
disp("Final Values:")
assumption = mod(assumption,360)
disp("Iterations taken: ")
disp(i)
calculated = double(subs(zeroTsix,symbols,assumption))

%plotting the skeleton
IRB_1410 = double(subs(IRB_1410,symbols,assumption))
x = linspace(-1500,1500,3001);
y = linspace(-1500,1500,3001);
[X,Y]=meshgrid(x,y);
HT= eye(4);
x_start = 0; y_start = 0; z_start = 0;
for i=1:6
    HT = HT * homoFromDH(IRB_1410(i,1),IRB_1410(i,2),IRB_1410(i,3),IRB_1410(i,4));
    plot3([x_start HT(13)],[y_start HT(14)],[z_start HT(15)],'LineWidth',3);grid on;hold on
    xlabel("x-axis");ylabel("y-axis");zlabel("z-axis");
    x_start = HT(13);y_start = HT(14);z_start = HT(15);
end
hold off
function mat = zRotHomo(ang)
    mat = [cosd(ang) -sind(ang) 0 0; sind(ang) cosd(ang) 0 0;0 0 1 0;0 0 0 1];
end

function mat = xRotHomo(ang)
    mat = [1 0 0 0;0 cosd(ang) -sind(ang) 0 ;0 sind(ang) cosd(ang) 0;0 0 0 1];
end

function homo = homoFromDH(a,alpha,d,theta)
    homo = transl(0,0,d)*zRotHomo(theta)*transl(a,0,0)*xRotHomo(alpha);
end