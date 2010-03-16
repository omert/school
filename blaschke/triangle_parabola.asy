texpreamble ("\usepackage{amsfonts}");

void
lineBarCircle(pair x1, pair x2, real R){
  real alpha = atan((x1.y-x2.y) / (x1.x-x2.x));
  draw(x1 + R * (cos(alpha), sin(alpha)) -- x2 - R * (cos(alpha), sin(alpha)));
}

void
triangleLineBackground(pair x0, pair x1, pair x2, int N){
  for (int i = 1; i < N; ++i){
    pair x1p = x0 + 1.0 * i / N * (x1 - x0);
    pair x2p = x0 + 1.0 * i / N * (x2 - x0);
    draw(x1p -- x2p);
  }
}

size (4cm);
real R = 0.025;
pair x0 = (0.0, 0.0);
pair x1 = (0.3, 1.0);
pair x3 = (1.0, 0.0);
pair x4 = x1 * 0.62;
pair x5 = x3 * 0.62;
pair x2 = (x4 + x5) / 2;


draw(x0 -- x1);
draw(x0 -- x3);
draw(x1 -- x3);
draw(x4 -- x2);
draw(x2 -- x5);
draw(x2 -- x1);
draw(x2 -- x3);
triangleLineBackground(x0, x1, x3, 30);
triangleLineBackground(x2, x4, x1, 10);
triangleLineBackground(x2, x5, x3, 10);

draw(x3{left} ..{x1}x1);

filldraw(circle(x2, R), white);
filldraw(circle(x0, R), white);
filldraw(circle(x1, R), white);
filldraw(circle(x3, R), white);
filldraw(circle(x4, R), white);
filldraw(circle(x5, R), white);

label("3", x3, E);
label("1", x1, E);
label("2", x2 + (0.06, 0.06));
fill(ellipse(x2 + (0.06, 0.06), 0.035, 0.05), white);
