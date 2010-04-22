import patterns;
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

real margin = 0.015 * Sx;
filldraw(currentpicture,(margin, K05 * Sy)--(K100 * Sx - margin, K05 * Sy)..(K141 * Sx - margin, K * Sy)..(K240y * Sx - margin, K240 * Sy)--(margin, K240 * Sy)--cycle, pattern("hatch"));

draw((0,0) -- (K20000 * Sx, K20000 * Sy));
draw((K240 * Sx - 0.05 * Sx, K240 * Sy) -- (K20000 * Sx - 0.05 * Sx, K20000 * Sy));
label("$A(w)=Kl(w)$",(K100 * Sx - 0.15 * Sx, K100 * Sy - 0.04 * Sy));
label("$A(w)={l(w)^2\over 20,000}$",(K20000 * Sx + 0.07 * Sx, K20000 * Sy - 0.1 * Sy));
draw((0,0)..(K100 * Sx, K05 * Sy)..(K141 * Sx, K * Sy)..(K240y * Sx, K240 * Sy)..(K20000 * Sx, K20000 * Sy));

pair w = (K20000 * Sx - 0.07 * Sx, K20000 * Sy - 0.05 * Sy);
label("$w$", w, SW);
draw(circle(w, 0.001* Sx));
label("FORBIDDEN ZONE", (K * Sx, K141* Sy));

draw((K100 * Sx, 0) -- (K100 * Sx, K100 * Sy));
