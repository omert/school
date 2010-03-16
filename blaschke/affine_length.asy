texpreamble ("\usepackage{amsfonts}");

void
triangle(pair x1, pair x2, real R, real r){
  pair x3 = (x1 + x2) / 2 * (1 + 1.0/2.0/1.732050808);
  draw(x1 -- x2);
  draw(x2 -- x3);
  draw(x3 -- x1);
  int N = 14;
  pair xv = x3 - (x1 + x2) / 2;
  for (int i = 1; i < N; ++i){
    pair xi1 = x1 + (x2 - x1) / N * i;
    pair xi2 = xi1 + xv * (1- abs(i - N/2) / N * 2);
    draw(xi1 -- xi2);
  }
  filldraw(circle(x1, r), white);
  filldraw(circle(x2, r), white);
  filldraw(circle(x3, r), white);
}

size (4cm);
real r = 0.05;
real R = 1;
real alpha0 = 95  * pi / 180;
real alpha1 = 155  * pi / 180;
real alpha2 = 215  * pi / 180;
real alpha3 = 275  * pi / 180;
real alphan = 335  * pi / 180;
pair x0 = R * (cos(alpha0), sin(alpha0));
pair x1 = R * (cos(alpha1), sin(alpha1));
pair x2 = R * (cos(alpha2), sin(alpha2));
pair x3 = R * (cos(alpha3), sin(alpha3));
pair xn = R * (cos(alphan), sin(alphan));

draw(arc((0.0, 0.0), R, 95, 360 / 6 * 4 + 95));

triangle(x0, x1, R, r);
triangle(x1, x2, R, r);
triangle(x2, x3, R, r);
triangle(x3, xn, R, r);
label("$\mathfrak{x}_0$", x0, E);
label("$\mathfrak{x}_1$", x1, SE);
label("$\mathfrak{x}_2$", x2, NE);
label("$\mathfrak{x}^0$", xn, N);
