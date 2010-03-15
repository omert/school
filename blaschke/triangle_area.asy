texpreamble ("\usepackage{amsfonts}");

void
lineBarCircle(real x1, real y1, real x2, real y2, real R, real theta){
  real x1p = x1 + y1 * cos(theta);
  real y1p = y1 * sin(theta);
  real x2p = x2 + y2 * cos(theta);
  real y2p = y2 * sin(theta);
  real alpha = atan((y2p-y1p) / (x2p-x1p));
  draw((x1p, y1p) + R * (cos(alpha), sin(alpha)) -- (x2p, y2p) - R * (cos(alpha), sin(alpha)));
}
void
pointWithLine(real x1, real x2, real R, real theta){
  real x1p = x1 + x2 * cos(theta);
  real x2p = x2 * sin(theta);
  draw (circle((x1p, x2p), R));
  //  label(L, (x1p+0.1, x2p));
  lineBarCircle(0.0, 0.0, x1, x2, R, theta);
  //  draw(R * (cos(alpha), sin(alpha)) -- (x1p, x2p) - R * (cos(alpha), sin(alpha)));
}

size (3cm);
real R = 0.04;
real x2len = 0.8;
draw (circle((0, 0), R));
real theta = pi/180*90;
pointWithLine(1.0, 0.0, R, theta);
draw((R,0.0)--(0.4,0.0),Arrow);
label("$\mathfrak{x}$",(0, -0.13));
label("$\mathfrak{x}'$",(0.4, -0.13));
label("$\mathfrak{z}$",(1, -0.13));

pointWithLine(0.1, 0.8, R, theta);
label("$\mathfrak{y}$",(0.0, 0.75));
lineBarCircle(0.1, 0.8, 1.0, 0.0, R, theta);
real alpha = atan(-0.8 / 0.9);
draw((0.1, 0.8) - R * (cos(alpha), sin(alpha)) -- (0.1, 0.8) - 10 * R * (cos(alpha), sin(alpha)), Arrow);
label("$\mathfrak{y}'$",(-0.25, 0.92));
