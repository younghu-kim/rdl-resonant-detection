\\ [사이클 #264] GL(3) A-gap 상관 — PARI/GP native
\\ δ-추출: Λ'/Λ(center+δ+it₀) → A = S₁² + 2H₁
\\ Sym²(11a1), N=121, d=3, center=1.5

default(realprecision, 38);
allocatemem(512000000);

\\ L-함수 초기화
E11 = ellinit([0,-1,1,-10,-10]);
L11 = lfunsympow(E11, 2);

\\ center 확인
k_center = L11[4]; \\ should be 3 (weight), center = 3/2 = 1.5
CENTER = k_center / 2;
print("center = ", CENTER);

\\ lfuninit — 넓은 범위
print("lfuninit 시작...");
L11i = lfuninit(L11, [0, 37]);
print("lfuninit OK");

\\ 영점 수집
zv = lfunzeros(L11i, 35);
nz = #zv;
print("영점: ", nz, "개");

\\ δ값
d1 = 0.005;
d2 = 0.01;
d3 = 0.02;

\\ 결과 파일 출력
outfile = "/tmp/gl3_c264_raw.csv";
write(outfile, "t,A_d1,A_d2,A_d3,S1_d1,H1_d1");

n_ok = 0;
n_fail = 0;

for(i = 1, nz,
  t0 = zv[i];
  if(t0 < 3.0, next);

  \\ δ₁
  ss1 = CENTER + d1 + t0*I;
  Lv1 = lfunlambda(L11i, ss1);
  if(abs(Lv1) < 1e-30,
    n_fail++; next
  );
  Ld1 = lfunlambda(L11i, ss1, 1);
  LpL1 = Ld1/Lv1;
  S1_1 = -imag(LpL1);
  H1_1 = (real(LpL1) - 1/d1) / d1;
  A1 = S1_1^2 + 2*H1_1;

  \\ δ₂
  ss2 = CENTER + d2 + t0*I;
  Lv2 = lfunlambda(L11i, ss2);
  if(abs(Lv2) < 1e-30,
    n_fail++; next
  );
  Ld2 = lfunlambda(L11i, ss2, 1);
  LpL2 = Ld2/Lv2;
  S1_2 = -imag(LpL2);
  H1_2 = (real(LpL2) - 1/d2) / d2;
  A2 = S1_2^2 + 2*H1_2;

  \\ δ₃
  ss3 = CENTER + d3 + t0*I;
  Lv3 = lfunlambda(L11i, ss3);
  if(abs(Lv3) < 1e-30,
    n_fail++; next
  );
  Ld3 = lfunlambda(L11i, ss3, 1);
  LpL3 = Ld3/Lv3;
  S1_3 = -imag(LpL3);
  H1_3 = (real(LpL3) - 1/d3) / d3;
  A3 = S1_3^2 + 2*H1_3;

  write(outfile, t0, ",", A1, ",", A2, ",", A3, ",", S1_1, ",", H1_1);
  n_ok++;

  if(i % 20 == 0, print("  ", i, "/", nz, " (", n_ok, " ok, ", n_fail, " fail)"));
);

print("Sym²(11a1) 완료: ", n_ok, " ok, ", n_fail, " fail");

\\ Sym²(37a1)
print("\n=== Sym²(37a1) ===");
E37 = ellinit([0,0,1,-1,0]);
L37 = lfunsympow(E37, 2);
L37i = lfuninit(L37, [0, 37]);
zv37 = lfunzeros(L37i, 35);
nz37 = #zv37;
print("영점: ", nz37, "개");

outfile37 = "/tmp/gl3_c264_37a1.csv";
write(outfile37, "t,A_d1,A_d2,A_d3,S1_d1,H1_d1");

n_ok37 = 0;
n_fail37 = 0;

for(i = 1, nz37,
  t0 = zv37[i];
  if(t0 < 3.0, next);

  ss1 = CENTER + d1 + t0*I;
  Lv1 = lfunlambda(L37i, ss1);
  if(abs(Lv1) < 1e-30, n_fail37++; next);
  Ld1 = lfunlambda(L37i, ss1, 1);
  LpL1 = Ld1/Lv1;
  S1_1 = -imag(LpL1);
  H1_1 = (real(LpL1) - 1/d1) / d1;
  A1 = S1_1^2 + 2*H1_1;

  ss2 = CENTER + d2 + t0*I;
  Lv2 = lfunlambda(L37i, ss2);
  if(abs(Lv2) < 1e-30, n_fail37++; next);
  Ld2 = lfunlambda(L37i, ss2, 1);
  LpL2 = Ld2/Lv2;
  A2 = (-imag(LpL2))^2 + 2*(real(LpL2) - 1/d2)/d2;

  ss3 = CENTER + d3 + t0*I;
  Lv3 = lfunlambda(L37i, ss3);
  if(abs(Lv3) < 1e-30, n_fail37++; next);
  Ld3 = lfunlambda(L37i, ss3, 1);
  LpL3 = Ld3/Lv3;
  A3 = (-imag(LpL3))^2 + 2*(real(LpL3) - 1/d3)/d3;

  write(outfile37, t0, ",", A1, ",", A2, ",", A3, ",", S1_1, ",", H1_1);
  n_ok37++;
  if(i % 20 == 0, print("  ", i, "/", nz37, " (", n_ok37, " ok)"));
);

print("Sym²(37a1) 완료: ", n_ok37, " ok, ", n_fail37, " fail");
print("\nDONE");
quit;
