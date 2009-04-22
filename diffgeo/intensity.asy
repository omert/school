import graph;
real R;
R = 100;
draw(arc((0, 0), R, 90, 270));
dot((-0.5 * R, 0));
label("$C$", (-0.5 * R, 0), NE);
draw((-R * 1.2, 0)--(R * 0.5, 0), EndArrow);
draw((0, -R * 1.2)--(0, R * 1.2), EndArrow);

real phi = pi * 5/6;
real rcphi = R * cos(phi);
real rsphi = R * sin(phi);
real yint = R * (sin(phi) - (0.5 + cos(phi)) * sin(2 * phi) / cos(2 * phi));
draw((R * 0.5, rsphi)--(rcphi, rsphi), MidArrow);
draw((rcphi, rsphi)--(-0.5 * R, yint),MidArrow);
draw((0, 0)--(rcphi, rsphi));
draw(arc((0, 0), R / 10, 0, phi/pi*180));
label("$\phi$", (R / 13, R / 7));
dot((-0.5 * R, yint));
label("$Y(y)$", (-0.5 * R, yint), E);

draw((-0.5 * R, 0){left}..(-0.61 * R, 0.07 * R));
