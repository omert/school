/*
texpreamble ("\usepackage{amsfonts}");

void
segment(pair x1, pair x2, real R, real r){
  draw(x1 -- x2);
  pair xb = x2 - x1;
  xb = (xb.y, -xb.x);
  int N = 14;
  
  for (int i = 1; i < N; ++i){
    pair xa = x1 + (x2 - x1) / N * i;
    real a = xb.x * xb.x + xb.y * xb.y;
    real b = 2 * (xa.x * xb.x + xa.y * xb.y);
    real c = xa.x * xa.x + xa.y * xa.y - R * R;
    real g = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    pair xc = xa + g * xb;
    draw(xa -- xc);
  }
  
  filldraw(circle(x1, r), white);
  filldraw(circle(x2, r), white);
}

size (4cm);
real r = 0.05;
real R = 1;
real alpha1 = 60  * pi / 180;
real alpha2 = 120  * pi / 180;
real alpha3 = 180  * pi / 180;
real alpha4 = 240  * pi / 180;
real alpha5 = 300  * pi / 180;
pair x1 = R * (cos(alpha1), sin(alpha1));
pair x2 = R * (cos(alpha2), sin(alpha2));
pair x3 = R * (cos(alpha3), sin(alpha3));
pair x4 = R * (cos(alpha4), sin(alpha4));
pair x5 = R * (cos(alpha5), sin(alpha5));

draw(arc((0.0, 0.0), R, 60, 300));

segment(x1, x2, R, r);
segment(x2, x3, R, r);
segment(x3, x4, R, r);
segment(x4, x5, R, r);

label("$\mathfrak{x}_1$", x1, E);
label("$\mathfrak{x}_2$", x2, NW);
label("$\mathfrak{x}_3$", x3, W);
label("$\mathfrak{x}_4$", x4, SW);
*/
import patterns;
//add("hatch",hatch(10mm, red));
add("hatch",checker(opacity(0.1)+red));

void
addLable(real pos, string decalx, string decaly, real Sx, real Sy){
  label(decalx, (pos * Sx, 0), S);
  label(decaly, (0, pos * Sy), W);
  draw((-0.005 * Sx, pos * Sy) -- (0.005 * Sx, pos * Sy));
  draw((pos * Sx, -0.005 * Sx) -- (pos * Sx, 0.005 * Sx));
}

real Sx = 400;
real Sy = 600;
draw((0,0) -- (0,Sy), Arrow);
label("$A(w)$", (0,Sy), E);
draw((0,0) -- (Sx,0), Arrow);
label("$l(w)$", (Sx,0), S);

real K05 = 0.15;
real K = 0.3;
real K100 = 0.45;
real K141 = 0.6;
real K240 = 0.75;
real K240y = 0.86;
real K20000 = 0.9;
addLable(K05, "${1\over 2}K$", "${1\over 2}K^2$", Sx, Sy);
addLable(K, "$K$", "$K^2$", Sx, Sy);
addLable(K100, "$100K$", "$100K^2$", Sx, Sy);
addLable(K141, "$141K$", "$141K^2$", Sx, Sy);
addLable(K240, "$240K$", "$240K^2$", Sx, Sy);
addLable(K20000, "$20,000K$", "$20,000K^2$", Sx, Sy);

real margin = 0.02 * Sx;
filldraw(currentpicture,(margin, K05 * Sy)--(K100 * Sx - margin, K05 * Sy)..(K141 * Sx - margin, K * Sy)..(K240y * Sx - margin, K240 * Sy)--(margin, K240 * Sy)--cycle, pattern("hatch"));

draw((0,0) -- (K20000 * Sx, K20000 * Sy));
label("$A(w)=l(w)$",(K20000 * Sx - 0.15 * Sx, K20000 * Sy - 0.04 * Sy));
label("$A(w)={l(w)^2\over 20,000}$",(K20000 * Sx + 0.1 * Sx, K20000 * Sy - 0.15 * Sy));
draw((0,0)..(K100 * Sx, K05 * Sy)..(K141 * Sx, K * Sy)..(K240y * Sx, K240 * Sy)..(K20000 * Sx, K20000 * Sy));

label("$w$", (K240 * Sx - 0.05 * Sx, K240 * Sy + 0.05 * Sy), SE);
draw(circle((K240 * Sx - 0.05 * Sx, K240 * Sy + 0.05 * Sy), 0.001* Sx));
label("FORBIDDEN ZONE", (K * Sx, K141* Sy));
