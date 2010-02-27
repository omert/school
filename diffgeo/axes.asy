texpreamble ("\usepackage{amsfonts}");
size (3cm);
real R = 0.04;
real x1 = 0.6;
real x2 = 0.5;
real x2len = 0.8;
draw (circle((0, 0), R));
draw ((R, 0) -- (x1 - R, 0));
draw (circle((x1, 0), R));
draw ((x1 + R, 0) -- (1,0), Arrow);
real theta = pi/180*80;
draw (R * (cos(theta), sin(theta)) -- x2len * (cos(theta), sin(theta)), Arrow);

draw (circle((x1 + x2 * cos(theta), x2 * sin(theta)), R));
draw ((x1 + R * cos(theta), R * sin(theta)) -- (x1 + (x2 - R) * cos(theta), (x2 - R) * sin(theta)));
label("$x_1$", (1, -0.1));
label("$x_2$", (x2len * cos(theta)-0.11, x2len * sin(theta)));
label("$\mathfrak{x}$",(x1+0.2, x2));
label("$\mathfrak{o}$",(0, -0.1));
