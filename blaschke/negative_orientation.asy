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
pointWithLine(real x1, real x2, real R, real theta, string L){
  real x1p = x1 + x2 * cos(theta);
  real x2p = x2 * sin(theta);
  draw (circle((x1p, x2p), R));
  label(L, (x1p+0.1, x2p));
  lineBarCircle(0.0, 0.0, x1, x2, R, theta);
  //  draw(R * (cos(alpha), sin(alpha)) -- (x1p, x2p) - R * (cos(alpha), sin(alpha)));
}

size (3cm);
real R = 0.04;
real x2len = 0.8;
draw (circle((0, 0), R));
draw ((R, 0) -- (1,0), Arrow);
real theta = pi/180*80;
draw (R * (cos(theta), sin(theta)) -- x2len * (cos(theta), sin(theta)), Arrow);

label("$x_1$", (1, -0.1));
label("$x_2$", (x2len * cos(theta)-0.11, x2len * sin(theta)));
label("$\mathfrak{o}$",(0, -0.1));

pointWithLine(0.6, 0.2, R, theta, "$\mathfrak{y}$");
pointWithLine(0.3, 0.6, R, theta, "$\mathfrak{x}$");
lineBarCircle(0.3, 0.6, 0.6, 0.2, R, theta);

