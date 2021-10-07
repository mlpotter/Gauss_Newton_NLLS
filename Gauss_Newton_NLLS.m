close all; clear all; clc;
N =7; r = 5;

theta = ((2 * pi / N) * (1:(N))).';
centerX = 0 ; centerY = 0;
x = centerX + r .* cos(theta);
y = centerY + r .* sin(theta);
beacons = [x y] ;
beacons = beacons + normrnd(0,.001,size(beacons)); 

theta = 2*pi*rand();
R = r/1.3*sqrt(rand());
centerX = 0 ; centerY = 0;
x = centerX + R  .* cos(theta);
y = centerY + R  .* sin(theta);
position = [x y];


% position = randi([-5,5],1,2);
distances = sqrt(sum((beacons-position).^2,2));
distances = distances + normrnd(0,.1,size(distances));
xk = randi([-15,15],1,2); %rk = distances - sqrt(sum((beacons-xk).^2,2));
%%
plot(beacons(:,1),beacons(:,2),'r*')
hold on; 
plot(position(1),position(2),'bo','markersize',5)
text(position(1),position(2),'True Position')
plot(xk(1),xk(2),'mx','markersize',5)
text(xk(1),xk(2),'Initial Position Estimation')

%%
mu =8;
for i = 1:25
    rk = residual(xk,beacons,distances);
    Ak = derivative(xk,beacons);
    bk = bias(Ak,xk,rk);
    xk = ( inv(Ak.'*Ak + mu*eye(size(Ak,2)) ) * (Ak' * bk + mu*xk.') )';
    plot(xk(1),xk(2),'gx')
    text(xk(1),xk(2),string(i))
    
%     if mod(i,1) == 0
%         position = [.2,.2] + position
%         distances = sqrt(sum((beacons-position).^2,2));
%         distances = distances + normrnd(0,.1,size(distances));
%         plot(position(1),position(2),'bo','markersize',5)
%     end
end
%%
%%
function A = derivative(xk,beacons)
   A = -(xk-beacons) ./  ( sqrt(sum((xk-beacons).^2,2)) + eps );
end

function b = bias(A,xk,rk)
   b = A*xk.' - rk;
end

function r = residual(xk,beacons,distances)
    r = distances - sqrt(sum((xk-beacons).^2,2));
end

