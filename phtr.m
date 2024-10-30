G = 6.67e-8;
M = 2.864e33;
c = 2.998e10;
rs = 2*G*M/(c^2);
function ua = ua(phi,u)
    ua = [u(2);3*(u(1).^2) - u(1)];
end
rm= 1e8;
rpuls = 1e6;
minlam = acos(sqrt(rpuls/rm));
phi = minlam;
alpha = 0.01;
u0 = 0.5*rs/rpuls;
uv0 = -u0/tan(alpha-phi);
[phi,l] = ode45(@ua,[minlam,-3*pi],[u0;uv0]);
r = 0.5*rs./(l(:,1));
res = [r.*cos(phi),r.*sin(phi)];
plot(res(:,1),res(:,2))
