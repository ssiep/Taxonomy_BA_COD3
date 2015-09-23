function [x] = focusSpericalToCartesian(rr, th, ph)


nu =   rr.*cos(th);
mu1 = -rr.*sin(th).*sin(ph);
mu2 =  rr.*sin(th).*cos(ph);

x = [nu;mu1;mu2];

end