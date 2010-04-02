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

