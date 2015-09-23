function out = Cod3FocusTimeFor_SphericalParam
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,rr,th,ph,bb)
nu=rr*cos(th);
mu1=-rr*sin(th)*sin(ph);
mu2=rr*sin(th)*cos(ph);
dydt=[kmrgd(2);
-kmrgd(1)^3+mu2*kmrgd(1)+mu1+kmrgd(2)*(nu+bb*kmrgd(1)+kmrgd(1)^2);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Cod3FocusTimeFor_SphericalParam);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,rr,th,ph,bb)
jac=[ 0 , 1 ; kmrgd(2)*(2*kmrgd(1) + bb) - 3*kmrgd(1)^2 + rr*cos(ph)*sin(th) , kmrgd(1)*bb + rr*cos(th) + kmrgd(1)^2 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,rr,th,ph,bb)
jacp=[ 0 , 0 , 0 , 0 ; kmrgd(2)*cos(th) - sin(ph)*sin(th) + kmrgd(1)*cos(ph)*sin(th) , kmrgd(1)*rr*cos(ph)*cos(th) - rr*cos(th)*sin(ph) - kmrgd(2)*rr*sin(th) , - rr*cos(ph)*sin(th) - kmrgd(1)*rr*sin(ph)*sin(th) , kmrgd(1)*kmrgd(2) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,rr,th,ph,bb)
hess1=[ 0 , 0 ; 2*kmrgd(2) - 6*kmrgd(1) , 2*kmrgd(1) + bb ];
hess2=[ 0 , 0 ; 2*kmrgd(1) + bb , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,rr,th,ph,bb)
hessp1=[ 0 , 0 ; cos(ph)*sin(th) , cos(th) ];
hessp2=[ 0 , 0 ; rr*cos(ph)*cos(th) , -rr*sin(th) ];
hessp3=[ 0 , 0 ; -rr*sin(ph)*sin(th) , 0 ];
hessp4=[ 0 , 0 ; kmrgd(2) , kmrgd(1) ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,rr,th,ph,bb)
tens31=[ 0 , 0 ; -6 , 2 ];
tens32=[ 0 , 0 ; 2 , 0 ];
tens33=[ 0 , 0 ; 2 , 0 ];
tens34=[ 0 , 0 ; 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,rr,th,ph,bb)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,rr,th,ph,bb)
