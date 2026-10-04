/* Calib — classical (non-learned) camera calibration math.
   Model: OpenCV pinhole + radial-tangential (k1,k2,p1,p2,k3).
   Pixel convention: pixel centres at integer coordinates, (0,0) = centre of top-left pixel. */
(function () {
'use strict';
const C = {};
const zeros = (n, m) => Array.from({ length: n }, () => new Array(m).fill(0));
const eye = n => { const a = zeros(n, n); for (let i = 0; i < n; i++) a[i][i] = 1; return a; };
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const nrm = a => Math.hypot(a[0], a[1], a[2]);
const unit = a => { const n = nrm(a); return [a[0] / n, a[1] / n, a[2] / n]; };
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const neg = a => [-a[0], -a[1], -a[2]];
const mm = (A, B) => A.map(r => B[0].map((_, j) => r.reduce((s, v, k) => s + v * B[k][j], 0)));
const mv = (A, v) => A.map(r => r.reduce((s, x, k) => s + x * v[k], 0));
const T = A => A[0].map((_, j) => A.map(r => r[j]));
const det3 = M => M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
const sumsq = r => { let s = 0; for (const v of r) s += v * v; return s; };
Object.assign(C, { dot, cross, nrm, unit, sub, mm, mv, T, det3 });

// ---------- linear algebra ----------
C.jacobiEig = function (A) {
  const n = A.length, a = A.map(r => r.slice()), v = eye(n);
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0, diag = 0;
    for (let i = 0; i < n; i++) { diag += a[i][i] * a[i][i]; for (let j = i + 1; j < n; j++) off += a[i][j] * a[i][j]; }
    if (off <= 1e-32 * diag || off === 0) break;
    for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
      const apq = a[p][q]; if (Math.abs(apq) < 1e-300) continue;
      const th = (a[q][q] - a[p][p]) / (2 * apq);
      const t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = c * x - s * y; a[k][q] = s * x + c * y; }
      for (let k = 0; k < n; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = c * x - s * y; a[q][k] = s * x + c * y; }
      for (let k = 0; k < n; k++) { const x = v[k][p], y = v[k][q]; v[k][p] = c * x - s * y; v[k][q] = s * x + c * y; }
    }
  }
  return { vals: a.map((r, i) => r[i]), vecs: v };
};
C.minEigvec = function (A) { const e = C.jacobiEig(A); let k = 0; for (let i = 1; i < e.vals.length; i++) if (e.vals[i] < e.vals[k]) k = i; return e.vecs.map(r => r[k]); };
C.solve = function (A, b) {
  const n = A.length, M = A.map((r, i) => r.concat([b[i]]));
  for (let c = 0; c < n; c++) {
    let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    if (!(Math.abs(M[p][c]) > 1e-300)) return null;
    [M[c], M[p]] = [M[p], M[c]];
    for (let r = c + 1; r < n; r++) { const f = M[r][c] / M[c][c]; if (f) for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k]; }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) { let s = M[i][n]; for (let k = i + 1; k < n; k++) s -= M[i][k] * x[k]; x[i] = s / M[i][i]; }
  return x.every(isFinite) ? x : null;
};
C.inv = function (A) {
  const n = A.length, I = eye(n), M = A.map((r, i) => r.concat(I[i]));
  for (let c = 0; c < n; c++) {
    let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    if (!(Math.abs(M[p][c]) > 1e-300)) return null;
    [M[c], M[p]] = [M[p], M[c]];
    const d = M[c][c]; for (let k = 0; k < 2 * n; k++) M[c][k] /= d;
    for (let r = 0; r < n; r++) { if (r === c) continue; const f = M[r][c]; if (f) for (let k = 0; k < 2 * n; k++) M[r][k] -= f * M[c][k]; }
  }
  return M.map(r => r.slice(n));
};

// ---------- rotations ----------
C.rodrigues = function (r) {
  const th = Math.hypot(r[0], r[1], r[2]);
  if (th < 1e-10) { const K = [[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]], K2 = mm(K, K); return eye(3).map((row, i) => row.map((v, j) => v + K[i][j] + 0.5 * K2[i][j])); }
  const k = [r[0] / th, r[1] / th, r[2] / th], c = Math.cos(th), s = Math.sin(th), v = 1 - c;
  return [[c + k[0] * k[0] * v, k[0] * k[1] * v - k[2] * s, k[0] * k[2] * v + k[1] * s],
          [k[1] * k[0] * v + k[2] * s, c + k[1] * k[1] * v, k[1] * k[2] * v - k[0] * s],
          [k[2] * k[0] * v - k[1] * s, k[2] * k[1] * v + k[0] * s, c + k[2] * k[2] * v]];
};
C.rotToRvec = function (R) {
  const tr = R[0][0] + R[1][1] + R[2][2]; let w, x, y, z, S;
  if (tr > 0) { S = Math.sqrt(tr + 1) * 2; w = 0.25 * S; x = (R[2][1] - R[1][2]) / S; y = (R[0][2] - R[2][0]) / S; z = (R[1][0] - R[0][1]) / S; }
  else if (R[0][0] > R[1][1] && R[0][0] > R[2][2]) { S = Math.sqrt(1 + R[0][0] - R[1][1] - R[2][2]) * 2; w = (R[2][1] - R[1][2]) / S; x = 0.25 * S; y = (R[0][1] + R[1][0]) / S; z = (R[0][2] + R[2][0]) / S; }
  else if (R[1][1] > R[2][2]) { S = Math.sqrt(1 + R[1][1] - R[0][0] - R[2][2]) * 2; w = (R[0][2] - R[2][0]) / S; x = (R[0][1] + R[1][0]) / S; y = 0.25 * S; z = (R[1][2] + R[2][1]) / S; }
  else { S = Math.sqrt(1 + R[2][2] - R[0][0] - R[1][1]) * 2; w = (R[1][0] - R[0][1]) / S; x = (R[0][2] + R[2][0]) / S; y = (R[1][2] + R[2][1]) / S; z = 0.25 * S; }
  if (w < 0) { w = -w; x = -x; y = -y; z = -z; }
  const sn = Math.hypot(x, y, z); if (sn < 1e-15) return [0, 0, 0];
  const th = 2 * Math.atan2(sn, w); return [x / sn * th, y / sn * th, z / sn * th];
};
C.nearestRotation = function (M) {
  let X = M.map(r => r.slice());
  for (let i = 0; i < 60; i++) {
    const Xi = C.inv(X); if (!Xi) break; const XiT = T(Xi); let d = 0;
    X = X.map((r, a) => r.map((v, b) => { const nv = 0.5 * (v + XiT[a][b]); d += Math.abs(nv - v); return nv; }));
    if (d < 1e-15) break;
  }
  return X;
};
// R = Rz(c)·Ry(b)·Rx(a), degrees
C.rotToEuler = R => { const b = -Math.asin(Math.max(-1, Math.min(1, R[2][0]))); return [Math.atan2(R[2][1], R[2][2]), b, Math.atan2(R[1][0], R[0][0])].map(v => v * 180 / Math.PI); };
C.eulerToRot = (ad, bd, cd) => {
  const [a, b, c] = [ad, bd, cd].map(v => v * Math.PI / 180);
  const Rx = [[1, 0, 0], [0, Math.cos(a), -Math.sin(a)], [0, Math.sin(a), Math.cos(a)]];
  const Ry = [[Math.cos(b), 0, Math.sin(b)], [0, 1, 0], [-Math.sin(b), 0, Math.cos(b)]];
  const Rz = [[Math.cos(c), -Math.sin(c), 0], [Math.sin(c), Math.cos(c), 0], [0, 0, 1]];
  return mm(Rz, mm(Ry, Rx));
};
C.lookAt = (Cw, target) => {
  const z = unit(sub(target, Cw)); let x = cross(z, [0, 0, 1]); if (nrm(x) < 1e-9) x = [1, 0, 0]; x = unit(x);
  const y = cross(z, x), R = [x, y, z]; return { R, t: neg(mv(R, Cw)) };
};
C.camCenter = c => neg(mv(T(c.R), c.t));

// ---------- camera model ----------
C.defaultCam = (w, h) => {
  const f = Math.max(w, h), L = C.lookAt([4, -6, 4], [4, 4, 0]);
  return { w, h, fx: f, fy: f, cx: (w - 1) / 2, cy: (h - 1) / 2, k1: 0, k2: 0, p1: 0, p2: 0, k3: 0, R: L.R, t: L.t, hasPose: false };
};
C.cloneCam = c => Object.assign({}, c, { R: c.R.map(r => r.slice()), t: c.t.slice() });
C.distortN = (c, x, y) => {
  const r2 = x * x + y * y, rad = 1 + r2 * (c.k1 + r2 * (c.k2 + r2 * c.k3));
  return [x * rad + 2 * c.p1 * x * y + c.p2 * (r2 + 2 * x * x), y * rad + c.p1 * (r2 + 2 * y * y) + 2 * c.p2 * x * y];
};
C.undistortN = (c, xd, yd, x0, y0) => {
  let x = x0 === undefined ? xd : x0, y = y0 === undefined ? yd : y0;
  if (!c.k1 && !c.k2 && !c.k3 && !c.p1 && !c.p2) return [xd, yd];
  for (let i = 0; i < 40; i++) {
    const r2 = x * x + y * y, rad = 1 + r2 * (c.k1 + r2 * (c.k2 + r2 * c.k3)), Rr = c.k1 + r2 * (2 * c.k2 + 3 * c.k3 * r2);
    const ex = x * rad + 2 * c.p1 * x * y + c.p2 * (r2 + 2 * x * x) - xd, ey = y * rad + c.p1 * (r2 + 2 * y * y) + 2 * c.p2 * x * y - yd;
    const a = rad + 2 * x * x * Rr + 2 * c.p1 * y + 6 * c.p2 * x, b = 2 * x * y * Rr + 2 * c.p1 * x + 2 * c.p2 * y, d = rad + 2 * y * y * Rr + 6 * c.p1 * y + 2 * c.p2 * x;
    const det = a * d - b * b; if (!(Math.abs(det) > 1e-12)) return [NaN, NaN];
    const dx = (d * ex - b * ey) / det, dy = (a * ey - b * ex) / det; x -= dx; y -= dy;
    if (dx * dx + dy * dy < 1e-28) break;
  }
  return isFinite(x) && isFinite(y) ? [x, y] : [NaN, NaN];
};
C.undistortPixel = (c, u, v) => { const q = C.undistortN(c, (u - c.cx) / c.fx, (v - c.cy) / c.fy); return [c.fx * q[0] + c.cx, c.fy * q[1] + c.cy]; };
C.distortPixel = (c, u, v) => { const q = C.distortN(c, (u - c.cx) / c.fx, (v - c.cy) / c.fy); return [c.fx * q[0] + c.cx, c.fy * q[1] + c.cy]; };
C.project = (c, X, ideal) => {
  const R = c.R, t = c.t;
  const xc = R[0][0] * X[0] + R[0][1] * X[1] + R[0][2] * X[2] + t[0], yc = R[1][0] * X[0] + R[1][1] * X[1] + R[1][2] * X[2] + t[1], zc = R[2][0] * X[0] + R[2][1] * X[1] + R[2][2] * X[2] + t[2];
  if (!(zc > 1e-9)) return null;
  let x = xc / zc, y = yc / zc; if (!ideal) [x, y] = C.distortN(c, x, y);
  return [c.fx * x + c.cx, c.fy * y + c.cy];
};
// keep pixel-space distortion unchanged when f changes by factor s
C.rescaleDist = (c, s) => { c.k1 *= s * s; c.k2 *= s ** 4; c.k3 *= s ** 6; c.p1 *= s; c.p2 *= s; return c; };
C.radiusLimit = c => {
  let r = 0; for (const [u, v] of [[0, 0], [c.w - 1, 0], [0, c.h - 1], [c.w - 1, c.h - 1]]) { const q = C.undistortN(c, (u - c.cx) / c.fx, (v - c.cy) / c.fy); if (isFinite(q[0])) r = Math.max(r, Math.hypot(q[0], q[1])); }
  return (r || 2) * 1.3;
};

// ---------- Levenberg–Marquardt ----------
C.lm = function (fun, x0, opt) {
  opt = opt || {};
  const maxIt = opt.maxIter || 300, n = x0.length;
  let x = x0.slice(), r = fun(x), cost = sumsq(r), lambda = 1e-3, iters = 0; const m = r.length;
  const jac = (x) => {
    const J = Array.from({ length: m }, () => new Array(n));
    for (let j = 0; j < n; j++) {
      const h = (opt.rel || 1e-6) * Math.max(opt.scale ? opt.scale[j] : 1, Math.abs(x[j]));
      const xp = x.slice(), xm = x.slice(); xp[j] += h; xm[j] -= h;
      const rp = fun(xp), rm = fun(xm); for (let i = 0; i < m; i++) J[i][j] = (rp[i] - rm[i]) / (2 * h);
    }
    return J;
  };
  const normal = J => { const A = zeros(n, n), g = new Array(n).fill(0);
    for (let i = 0; i < m; i++) { const Ji = J[i]; for (let a = 0; a < n; a++) { g[a] += Ji[a] * r[i]; for (let b = a; b < n; b++) A[a][b] += Ji[a] * Ji[b]; } }
    for (let a = 0; a < n; a++) for (let b = 0; b < a; b++) A[a][b] = A[b][a]; return { A, g }; };
  let J = jac(x);
  for (; iters < maxIt; iters++) {
    const { A, g } = normal(J); let improved = false, conv = false;
    for (let tries = 0; tries < 30; tries++) {
      const Ad = A.map((row, i) => row.map((v, j) => i === j ? v + lambda * Math.max(v, 1e-12) : v));
      const d = C.solve(Ad, g.map(v => -v)); if (!d) { lambda *= 10; continue; }
      const xn = x.map((v, i) => v + d[i]), rn = fun(xn), cn = sumsq(rn);
      if (isFinite(cn) && cn < cost) {
        const old = cost; x = xn; r = rn; cost = cn; lambda = Math.max(lambda / 10, 1e-15); improved = true;
        conv = (old - cn) <= 1e-15 * old + 1e-30 || Math.hypot(...d) <= 1e-13 * (Math.hypot(...x) + 1e-13); break;
      }
      lambda *= 10; if (lambda > 1e16) break;
    }
    if (!improved || conv) break;
    J = jac(x);
  }
  J = jac(x); const { A } = normal(J), Ai = C.inv(A), dof = Math.max(1, m - n - (opt.extraDof || 0));
  const s2 = cost / dof, cov = Ai ? Ai.map(row => row.map(v => v * s2)) : null;
  return { x, r, cost, rms: Math.sqrt(cost / m), iters, cov, std: cov ? cov.map((row, i) => Math.sqrt(Math.max(0, row[i]))) : x.map(() => NaN), dof };
};

// ---------- normalisation ----------
function norm2T(P) { const n = P.length; let mx = 0, my = 0; P.forEach(p => { mx += p[0]; my += p[1]; }); mx /= n; my /= n; let d = 0; P.forEach(p => d += Math.hypot(p[0] - mx, p[1] - my)); d /= n; const s = d > 0 ? Math.SQRT2 / d : 1; return [[s, 0, -s * mx], [0, s, -s * my], [0, 0, 1]]; }
function norm3T(P) { const n = P.length; const m = [0, 0, 0]; P.forEach(p => { for (let k = 0; k < 3; k++) m[k] += p[k] / n; }); let d = 0; P.forEach(p => d += nrm(sub(p, m))); d /= n; const s = d > 0 ? Math.sqrt(3) / d : 1; return [[s, 0, 0, -s * m[0]], [0, s, 0, -s * m[1]], [0, 0, s, -s * m[2]], [0, 0, 0, 1]]; }

// ---------- homography (plane XY -> pixel) ----------
C.applyH = (H, p) => { const w = H[2][0] * p[0] + H[2][1] * p[1] + H[2][2]; if (!(w > 1e-12)) return null; return [(H[0][0] * p[0] + H[0][1] * p[1] + H[0][2]) / w, (H[1][0] * p[0] + H[1][1] * p[1] + H[1][2]) / w]; };
C.homography = function (src, dst) {
  const n = src.length; if (n < 4) return { ok: false, msg: 'Need ≥ 4 points.' };
  const Ts = norm2T(src), Td = norm2T(dst), A = zeros(9, 9);
  for (let i = 0; i < n; i++) {
    const s = mv(Ts, [src[i][0], src[i][1], 1]), d = mv(Td, [dst[i][0], dst[i][1], 1]), x = s[0], y = s[1], u = d[0], v = d[1];
    for (const row of [[-x, -y, -1, 0, 0, 0, u * x, u * y, u], [0, 0, 0, -x, -y, -1, v * x, v * y, v]])
      for (let a = 0; a < 9; a++) for (let b = 0; b < 9; b++) A[a][b] += row[a] * row[b];
  }
  const h = C.minEigvec(A), Hn = [[h[0], h[1], h[2]], [h[3], h[4], h[5]], [h[6], h[7], h[8]]], Tdi = C.inv(Td);
  let H = mm(Tdi, mm(Hn, Ts));
  let ws = 0; src.forEach(p => ws += H[2][0] * p[0] + H[2][1] * p[1] + H[2][2]);
  const sc = (ws < 0 ? -1 : 1) / Math.hypot(...H.flat()); H = H.map(r => r.map(v => v * sc));
  // geometric refinement (all 9 entries, gauge held by damping; renormalised after)
  const res = x => { const Hx = [x.slice(0, 3), x.slice(3, 6), x.slice(6, 9)], r = []; for (let i = 0; i < n; i++) { const p = C.applyH(Hx, src[i]); if (!p) { r.push(1e4, 1e4); continue; } r.push(p[0] - dst[i][0], p[1] - dst[i][1]); } return r; };
  const L = C.lm(res, H.flat(), { rel: 1e-7, extraDof: -1 });
  const nn = Math.hypot(...L.x); H = [L.x.slice(0, 3), L.x.slice(3, 6), L.x.slice(6, 9)].map(r => r.map(v => v / nn));
  const r = res(H.flat()), errs = []; for (let i = 0; i < n; i++) errs.push(Math.hypot(r[2 * i], r[2 * i + 1]));
  return { ok: true, H, rms: Math.sqrt(sumsq(r) / n), errs, dof: 2 * n - 8 };
};

// ---------- projection matrix DLT ----------
C.dlt = function (W, U) {
  const T2 = norm2T(U), T3 = norm3T(W), A = zeros(12, 12);
  for (let i = 0; i < W.length; i++) {
    const X = mv(T3, [W[i][0], W[i][1], W[i][2], 1]), x = mv(T2, [U[i][0], U[i][1], 1]), u = x[0], v = x[1];
    const r1 = [X[0], X[1], X[2], X[3], 0, 0, 0, 0, -u * X[0], -u * X[1], -u * X[2], -u * X[3]];
    const r2 = [0, 0, 0, 0, X[0], X[1], X[2], X[3], -v * X[0], -v * X[1], -v * X[2], -v * X[3]];
    for (const row of [r1, r2]) for (let a = 0; a < 12; a++) { if (!row[a]) continue; for (let b = 0; b < 12; b++) A[a][b] += row[a] * row[b]; }
  }
  const p = C.minEigvec(A), Pn = [p.slice(0, 4), p.slice(4, 8), p.slice(8, 12)];
  let P = mm(C.inv(T2), mm(Pn, T3));
  let s = 0; W.forEach(X => s += P[2][0] * X[0] + P[2][1] * X[1] + P[2][2] * X[2] + P[2][3]);
  if (s < 0) P = P.map(r => r.map(v => -v));
  return P;
};
C.decomposeP = function (P) {
  const sc = nrm([P[2][0], P[2][1], P[2][2]]); P = P.map(r => r.map(v => v / sc));
  const m1 = P[0].slice(0, 3), m2 = P[1].slice(0, 3), r3 = P[2].slice(0, 3);
  const cx = dot(m1, r3), cy = dot(m2, r3), u2 = sub(m2, r3.map(v => v * cy)), fy = nrm(u2), r2 = u2.map(v => v / fy);
  const sk = dot(m1, r2), u1 = sub(sub(m1, r3.map(v => v * cx)), r2.map(v => v * sk)), fx = nrm(u1), r1 = u1.map(v => v / fx);
  const R = [r1, r2, r3];
  if (det3(R) < 0) return { ok: false, msg: 'DLT gave a reflection: world axes appear left-handed (or the image is mirrored).' };
  const K = [[fx, sk, cx], [0, fy, cy], [0, 0, 1]], t = C.solve(K, [P[0][3], P[1][3], P[2][3]]);
  return { ok: true, K, R, t, skew: sk };
};

// ---------- plane fit ----------
C.planeFit = function (W) {
  const n = W.length, c = [0, 0, 0]; W.forEach(p => { for (let k = 0; k < 3; k++) c[k] += p[k] / n; });
  const S = zeros(3, 3); W.forEach(p => { const d = sub(p, c); for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) S[a][b] += d[a] * d[b] / n; });
  const e = C.jacobiEig(S), idx = [0, 1, 2].sort((a, b) => e.vals[b] - e.vals[a]), col = k => [e.vecs[0][k], e.vecs[1][k], e.vecs[2][k]];
  const l = idx.map(i => Math.max(0, e.vals[i])), e1 = unit(col(idx[0])), e2 = unit(col(idx[1])), e3 = cross(e1, e2);
  return { c, B: [e1, e2, e3], ratio: l[0] > 0 ? Math.sqrt(l[2] / l[0]) : 0, collinear: !(l[1] > 1e-12 * l[0]) };
};
C.zhangF = function (H, cx, cy) {
  let G = mm([[1, 0, -cx], [0, 1, -cy], [0, 0, 1]], H); const nn = Math.hypot(...G.flat()); G = G.map(r => r.map(v => v / nn));
  const h1 = [G[0][0], G[1][0], G[2][0]], h2 = [G[0][1], G[1][1], G[2][1]];
  const a1 = h1[0] * h2[0] + h1[1] * h2[1], b1 = h1[2] * h2[2], a2 = h1[0] ** 2 + h1[1] ** 2 - h2[0] ** 2 - h2[1] ** 2, b2 = h1[2] ** 2 - h2[2] ** 2;
  const w = -(a1 * b1 + a2 * b2) / (a1 * a1 + a2 * a2);
  if (!(w > 0) || !isFinite(w)) return { ok: false, msg: 'Plane homography gives no real focal length (plane too fronto-parallel or points inconsistent)' };
  return { ok: true, f: 1 / Math.sqrt(w) };
};
C.poseFromH = function (H, K) {
  const A = mm(C.inv(K), H), h1 = [A[0][0], A[1][0], A[2][0]], h2 = [A[0][1], A[1][1], A[2][1]], h3 = [A[0][2], A[1][2], A[2][2]];
  let lam = 2 / (nrm(h1) + nrm(h2)); if (h3[2] * lam < 0) lam = -lam;
  const r1 = h1.map(v => v * lam), r2 = h2.map(v => v * lam), t = h3.map(v => v * lam), r3 = cross(r1, r2);
  return { R: C.nearestRotation([[r1[0], r2[0], r3[0]], [r1[1], r2[1], r3[1]], [r1[2], r2[2], r3[2]]]), t };
};

// ---------- 2D-3D calibration ----------
// pts: [{img:[u,v], world:[X,Y,Z]}]; o: {estF, fixAspect, estPP, distModel, initCur}
C.DIST_MODELS = { keep: [], k1: ['k1'], k1k2: ['k1', 'k2'], k1k2p: ['k1', 'k2', 'p1', 'p2'], all: ['k1', 'k2', 'p1', 'p2', 'k3'] };
C.calibratePoints = function (pts, cam, o) {
  const warn = [], log = [], n = pts.length;
  if (n < 4) return { ok: false, msg: 'Need at least 4 points with world coordinates.' };
  const W = pts.map(p => p.world), I = pts.map(p => p.img), U = I.map(p => C.undistortPixel(cam, p[0], p[1]));
  if (U.some(p => !isFinite(p[0]))) return { ok: false, msg: 'Current distortion cannot be inverted at some points — reset or re-solve the lens.' };
  const pf = C.planeFit(W); if (pf.collinear) return { ok: false, msg: 'World points are collinear — the pose is undetermined.' };
  const c0 = C.cloneCam(cam); let estF = o.estF, estPP = o.estPP; const planar = pf.ratio < 0.02;
  if (estPP && planar) { warn.push('Principal point is not observable from a single planar target — held fixed.'); estPP = false; }
  if (o.initCur && cam.hasPose) log.push('Initialised from the current camera.');
  else if (planar) {
    const Q = W.map(X => { const d = sub(X, pf.c); return [dot(pf.B[0], d), dot(pf.B[1], d)]; });
    const h = C.homography(Q, U); if (!h.ok) return h;
    if (estF) { const z = C.zhangF(h.H, c0.cx, c0.cy); if (z.ok) { c0.fx = c0.fy = z.f; log.push('f initialised from plane homography: ' + z.f.toFixed(2) + ' px'); } else warn.push(z.msg + ' — started from current f.'); }
    const ps = C.poseFromH(h.H, [[c0.fx, 0, c0.cx], [0, c0.fy, c0.cy], [0, 0, 1]]);
    c0.R = mm(ps.R, pf.B); c0.t = sub(ps.t, mv(c0.R, pf.c));
    const tilt = Math.acos(Math.min(1, Math.abs(ps.R[2][2]))) * 180 / Math.PI;
    log.push('Pose initialised from homography (' + n + ' coplanar points, plane tilt ' + tilt.toFixed(1) + '°).');
    if (estF && tilt < 12) warn.push('Target plane is within ' + tilt.toFixed(1) + '° of fronto-parallel: focal length is weakly constrained.');
  } else if (n >= 6) {
    if (pf.ratio < 0.08) warn.push('Points are nearly coplanar — DLT is poorly conditioned; check σ values.');
    if (estF || estPP) {
      const d = C.decomposeP(C.dlt(W, U)); if (!d.ok) return d;
      c0.fx = d.K[0][0]; c0.fy = d.K[1][1]; if (o.fixAspect) c0.fx = c0.fy = Math.sqrt(c0.fx * c0.fy);
      if (!estF) { c0.fx = cam.fx; c0.fy = cam.fy; }
      if (estPP) { c0.cx = d.K[0][2]; c0.cy = d.K[1][2]; }
      c0.R = d.R; c0.t = d.t;
      log.push('DLT: fx ' + d.K[0][0].toFixed(1) + ', fy ' + d.K[1][1].toFixed(1) + ', c (' + d.K[0][2].toFixed(1) + ', ' + d.K[1][2].toFixed(1) + '), skew ' + d.skew.toFixed(2) + ' (dropped)');
    } else {
      const Xn = U.map(u => [(u[0] - cam.cx) / cam.fx, (u[1] - cam.cy) / cam.fy]), P = C.dlt(W, Xn), M = P.map(r => r.slice(0, 3)), dM = det3(M);
      if (!(dM > 0)) return { ok: false, msg: 'Calibrated DLT gave a reflection: world axes appear left-handed.' };
      const lam = Math.cbrt(dM); c0.R = C.nearestRotation(M.map(r => r.map(v => v / lam))); c0.t = [P[0][3] / lam, P[1][3] / lam, P[2][3] / lam];
      log.push('Pose initialised by calibrated DLT (K held).');
    }
  } else if (cam.hasPose) log.push('Fewer than 6 non-coplanar points: started from the current camera.');
  else return { ok: false, msg: 'Need ≥ 6 non-coplanar or ≥ 4 coplanar points — or set a starting pose in Verify.' };

  const names = [];
  if (estF) { if (o.fixAspect) names.push('f'); else names.push('fx', 'fy'); }
  if (estPP) names.push('cx', 'cy');
  const dm = C.DIST_MODELS[o.distModel] || []; names.push(...dm);
  const lock = dm.length === 0, f0 = Math.sqrt(cam.fx * cam.fy), base = { k1: cam.k1, k2: cam.k2, p1: cam.p1, p2: cam.p2, k3: cam.k3 };
  if (o.fixAspect && estF) c0.fy = c0.fx;
  const rv0 = C.rotToRvec(c0.R), np = names.length;
  const x0 = names.map(nm => nm === 'f' ? c0.fx : c0[nm]).concat(rv0, c0.t);
  const scale = names.map(nm => /^(f|fx|fy|cx|cy)$/.test(nm) ? 1 : 0.01).concat([0.01, 0.01, 0.01], c0.t.map(() => Math.max(1e-3, nrm(c0.t) * 1e-2)));
  const build = x => {
    const c = C.cloneCam(c0);
    names.forEach((nm, i) => { if (nm === 'f') { c.fx = x[i]; c.fy = x[i]; } else c[nm] = x[i]; });
    if (lock) { const s = Math.sqrt(c.fx * c.fy) / f0; Object.assign(c, base); C.rescaleDist(c, s); }
    c.R = C.rodrigues(x.slice(np, np + 3)); c.t = x.slice(np + 3, np + 6); return c;
  };
  if (2 * n <= np + 6) return { ok: false, msg: `Only ${2 * n} residuals for ${np + 6} parameters — add points or fix more parameters.` };
  const res = x => { const c = build(x), r = []; for (let i = 0; i < n; i++) { const p = C.project(c, W[i]); if (!p || !isFinite(p[0])) { r.push(1e4, 1e4); continue; } r.push(p[0] - I[i][0], p[1] - I[i][1]); } return r; };
  const L = C.lm(res, x0, { scale });
  const out = build(L.x); out.hasPose = true;
  const errs = []; for (let i = 0; i < n; i++) errs.push(Math.hypot(L.r[2 * i], L.r[2 * i + 1]));
  if (errs.some(e => e > 1e3)) warn.push('Some points project behind the camera — check world coordinates.');
  const std = {}; names.forEach((nm, i) => std[nm] = L.std[i]);
  std.rot_deg = Math.hypot(L.std[np], L.std[np + 1], L.std[np + 2]) * 180 / Math.PI; std.t = Math.hypot(L.std[np + 3], L.std[np + 4], L.std[np + 5]);
  if (2 * n < 2 * (np + 6)) warn.push(`Low redundancy: ${2 * n} residuals for ${np + 6} parameters. σ values are optimistic.`);
  return { ok: true, cam: out, rms: Math.sqrt(L.cost / n), errs, maxErr: Math.max(...errs), std, names, iters: L.iters, dof: L.dof, planar, warn, log };
};

// ---------- plumb-line distortion ----------
C.lineFitTLS = P => {
  const n = P.length; let mx = 0, my = 0; P.forEach(p => { mx += p[0] / n; my += p[1] / n; });
  let sxx = 0, sxy = 0, syy = 0; P.forEach(p => { const dx = p[0] - mx, dy = p[1] - my; sxx += dx * dx; sxy += dx * dy; syy += dy * dy; });
  const th = 0.5 * Math.atan2(2 * sxy, sxx - syy), nx = -Math.sin(th), ny = Math.cos(th);
  return { c: [mx, my], dir: [Math.cos(th), Math.sin(th)], n: [nx, ny], d: P.map(p => (p[0] - mx) * nx + (p[1] - my) * ny) };
};
C.polyStraightness = (pts, cam) => { if (pts.length < 3) return NaN; const u = pts.map(q => C.undistortPixel(cam, q[0], q[1])); if (u.some(q => !isFinite(q[0]))) return NaN; const d = C.lineFitTLS(u).d; return Math.sqrt(sumsq(d) / d.length); };
C.solveLens = function (polys, cam, model) {
  const L = polys.filter(p => p.length >= 3); if (!L.length) return { ok: false, msg: 'Trace at least one straight edge with ≥ 3 points.' };
  const names = { k1: ['k1'], k1k2: ['k1', 'k2'], k1k2k3: ['k1', 'k2', 'k3'] }[model] || ['k1', 'k2'];
  const dof = L.reduce((s, p) => s + p.length - 2, 0);
  if (dof <= names.length) return { ok: false, msg: `${dof} residual DOF for ${names.length} parameters — add points or lines.` };
  const make = x => { const c = Object.assign({}, cam, { k1: 0, k2: 0, k3: 0 }); names.forEach((nm, i) => c[nm] = x[i]); return c; };
  const res = x => { const c = make(x), r = []; for (const p of L) { const u = p.map(q => C.undistortPixel(c, q[0], q[1])); if (u.some(q => !isFinite(q[0]))) { p.forEach(() => r.push(1e3)); continue; } C.lineFitTLS(u).d.forEach(d => r.push(d)); } return r; };
  const zero = names.map(() => 0), r0 = res(zero), m = r0.length, before = Math.sqrt(sumsq(r0) / m);
  let best = null;
  for (const start of [zero, names.map(nm => cam[nm] || 0)]) { const R = C.lm(res, start, { scale: names.map(() => 0.01), extraDof: 2 * L.length }); if (!best || R.cost < best.cost) best = R; }
  const vals = {}, std = {}; names.forEach((nm, i) => { vals[nm] = best.x[i]; std[nm] = best.std[i]; });
  const c = make(best.x), warn = [];
  const rl = C.radiusLimit(c);
  // monotonicity of the radial map over the image (otherwise the model folds)
  for (let r = 0; r <= rl; r += rl / 200) { const r2 = r * r, d = 1 + 3 * c.k1 * r2 + 5 * c.k2 * r2 * r2 + 7 * c.k3 * r2 ** 3; if (d <= 0) { warn.push('Radial model folds inside the image — result is unreliable; add lines nearer the borders or use a simpler model.'); break; } }
  const sel = L.map(p => { const u = p.map(q => Math.hypot(q[0] - cam.cx, q[1] - cam.cy)); return Math.max(...u); });
  if (Math.max(...sel) < 0.35 * Math.hypot(cam.w, cam.h) / 2) warn.push('All lines are near the image centre where distortion is weak — add edges near the borders.');
  return { ok: true, vals, std, names, before, after: Math.sqrt(best.cost / m), perLine: L.map(p => C.polyStraightness(p, c)), warn, n: m };
};

// ---------- vanishing points ----------
C.fitVP = function (segs, s, c0) {
  const N = p => [(p[0] - c0[0]) * s, (p[1] - c0[1]) * s, 1], M = zeros(3, 3);
  for (const [a, b] of segs) { let l = cross(N(a), N(b)); const nn = Math.hypot(l[0], l[1]); l = l.map(v => v / nn); const w = Math.hypot(b[0] - a[0], b[1] - a[1]) * s; for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) M[i][j] += w * l[i] * l[j]; }
  let v = unit(C.minEigvec(M)); if (v[2] < 0) v = neg(v);
  const sph = x => [Math.sin(x[0]) * Math.cos(x[1]), Math.sin(x[0]) * Math.sin(x[1]), Math.cos(x[0])];
  const res = x => { const vv = sph(x), r = []; for (const [a, b] of segs) { const A = N(a), B = N(b), m = [(A[0] + B[0]) / 2, (A[1] + B[1]) / 2, 1]; let L = cross(m, vv); const nn = Math.hypot(L[0], L[1]); if (!(nn > 1e-15)) { r.push(0); continue; } r.push(dot(L, A) / nn / s); } return r; };
  const L = C.lm(res, [Math.acos(Math.max(-1, Math.min(1, v[2]))), Math.atan2(v[1], v[0])], { scale: [0.01, 0.01] });
  v = sph(L.x); if (v[2] < 0) v = neg(v);
  const finite = Math.abs(v[2]) > 1e-9, pix = finite ? [v[0] / v[2] / s + c0[0], v[1] / v[2] / s + c0[1]] : null;
  const ang = segs.map(([a, b]) => { const m = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2], d1 = [b[0] - a[0], b[1] - a[1]], d2 = finite ? [pix[0] - m[0], pix[1] - m[1]] : [v[0], v[1]]; const c = Math.abs(d1[0] * d2[0] + d1[1] * d2[1]) / (Math.hypot(...d1) * Math.hypot(...d2)); return Math.acos(Math.min(1, c)) * 180 / Math.PI; });
  return { v, pix, finite, rmsPx: L.rms, ang };
};
C.calibFromVPs = function (vps, axes, mode, pp) {
  const pairs = []; for (let i = 0; i < axes.length; i++) for (let j = i + 1; j < axes.length; j++) pairs.push([vps[axes[i]].v, vps[axes[j]].v]);
  let cx, cy, f2;
  if (mode === 'ortho') {
    if (axes.length < 3) return { ok: false, msg: 'Orthocentre mode needs lines on all three axes.' };
    const A = zeros(4, 4); for (const [a, b] of pairs) { const row = [a[0] * b[0] + a[1] * b[1], a[0] * b[2] + a[2] * b[0], a[1] * b[2] + a[2] * b[1], a[2] * b[2]]; for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) A[i][j] += row[i] * row[j]; }
    const w = C.minEigvec(A); if (Math.abs(w[0]) < 1e-12) return { ok: false, msg: 'Degenerate VP configuration.' };
    cx = -w[1] / w[0]; cy = -w[2] / w[0]; f2 = w[3] / w[0] - cx * cx - cy * cy;
  } else {
    cx = pp[0]; cy = pp[1]; let num = 0, den = 0;
    for (const [a, b] of pairs) { const A = a[2] * b[2], B = (a[0] - cx * a[2]) * (b[0] - cx * b[2]) + (a[1] - cy * a[2]) * (b[1] - cy * b[2]); num += A * B; den += A * A; }
    if (den < 1e-10) return { ok: false, msg: 'Vanishing points are at (or near) infinity: focal length is unobservable. Use a view with stronger perspective, or Points.' };
    f2 = -num / den;
  }
  if (!(f2 > 0)) return { ok: false, msg: 'Inconsistent VPs — no real focal length. Check that each axis group is parallel in 3D and the groups are mutually orthogonal.' };
  const f = Math.sqrt(f2), d = {};
  axes.forEach(k => { const v = vps[k].v; d[k] = unit([(v[0] - cx * v[2]) / f, (v[1] - cy * v[2]) / f, v[2]]); });
  const orth = []; for (let i = 0; i < axes.length; i++) for (let j = i + 1; j < axes.length; j++) orth.push(Math.abs(90 - Math.acos(Math.min(1, Math.abs(dot(d[axes[i]], d[axes[j]])))) * 180 / Math.PI));
  const r = Object.assign({}, d);
  if (r.z && r.z[1] > 0) r.z = neg(r.z);
  if (r.x && r.y) { if (r.y[2] < 0) r.y = neg(r.y); const zc = cross(r.x, r.y); if (r.z) { if (dot(zc, r.z) < 0) r.x = neg(r.x); } else { if (zc[1] > 0) r.x = neg(r.x); r.z = cross(r.x, r.y); } }
  else if (r.x && r.z) { if (r.x[2] < 0) r.x = neg(r.x); r.y = cross(r.z, r.x); }
  else if (r.y && r.z) { if (r.y[2] < 0) r.y = neg(r.y); r.x = cross(r.y, r.z); }
  const R = C.nearestRotation([[r.x[0], r.y[0], r.z[0]], [r.x[1], r.y[1], r.z[1]], [r.x[2], r.y[2], r.z[2]]]);
  return { ok: true, f, cx, cy, R, orth };
};
function rng(seed) { let a = seed >>> 0; return () => { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
C.solveVP = function (segsByAxis, cam, o) {
  const s = 1 / Math.max(cam.w, cam.h), c0 = [(cam.w - 1) / 2, (cam.h - 1) / 2];
  const und = {}; for (const k of ['x', 'y', 'z']) und[k] = (segsByAxis[k] || []).map(([a, b]) => [C.undistortPixel(cam, a[0], a[1]), C.undistortPixel(cam, b[0], b[1])]);
  const axes = ['x', 'y', 'z'].filter(k => und[k].length >= 2);
  if (axes.length < 2) return { ok: false, msg: 'Draw ≥ 2 lines on at least two axes.' };
  const ppC = o.ppMode === 'center' ? [0, 0] : [(cam.cx - c0[0]) * s, (cam.cy - c0[1]) * s];
  const run = U => { const vps = {}; axes.forEach(k => vps[k] = C.fitVP(U[k], s, c0)); const cal = C.calibFromVPs(vps, axes, o.ppMode, ppC); return { vps, cal }; };
  const main = run(und); if (!main.cal.ok) return { ok: false, msg: main.cal.msg, vps: main.vps };
  let sig2 = 0, cnt = 0; axes.forEach(k => { sig2 += main.vps[k].rmsPx ** 2 * und[k].length; cnt += und[k].length; });
  const sigma = Math.max(0.25, Math.sqrt(sig2 / cnt)), rnd = rng(12345), gauss = () => Math.sqrt(-2 * Math.log(1 - rnd())) * Math.cos(2 * Math.PI * rnd());
  const F = [], CX = [], CY = [];
  for (let i = 0; i < 120; i++) {
    const U = {}; axes.forEach(k => U[k] = und[k].map(([a, b]) => [[a[0] + sigma * gauss(), a[1] + sigma * gauss()], [b[0] + sigma * gauss(), b[1] + sigma * gauss()]]));
    const r = run(U); if (r.cal.ok) { F.push(r.cal.f / s); CX.push(r.cal.cx / s + c0[0]); CY.push(r.cal.cy / s + c0[1]); }
  }
  const sd = a => { if (a.length < 3) return NaN; const m = a.reduce((x, y) => x + y, 0) / a.length; return Math.sqrt(a.reduce((x, y) => x + (y - m) ** 2, 0) / (a.length - 1)); };
  let R = main.cal.R; const k = ((o.axisRot || 0) % 4 + 4) % 4; if (k) R = mm(R, C.eulerToRot(0, 0, 90 * k));
  const vpsPix = {}; axes.forEach(a => vpsPix[a] = { pix: main.vps[a].pix, dir: main.vps[a].v, finite: main.vps[a].finite, rmsPx: main.vps[a].rmsPx, ang: main.vps[a].ang });
  const warn = []; if (F.length < 100) warn.push('Solution is unstable under ' + sigma.toFixed(2) + ' px endpoint noise (' + (120 - F.length) + '/120 trials failed).');
  const fS = sd(F); if (fS / (main.cal.f / s) > 0.05) warn.push('Focal length uncertainty above 5% — use longer lines that converge more strongly.');
  if (main.cal.orth.length && Math.max(...main.cal.orth) > 1) warn.push('Recovered axes deviate ' + Math.max(...main.cal.orth).toFixed(2) + '° from orthogonal before projection onto SO(3).');
  return { ok: true, f: main.cal.f / s, cx: main.cal.cx / s + c0[0], cy: main.cal.cy / s + c0[1], R, vps: vpsPix, axes, sigma, std: { f: fS, cx: o.ppMode === 'ortho' ? sd(CX) : 0, cy: o.ppMode === 'ortho' ? sd(CY) : 0 }, orth: main.cal.orth, warn };
};

// ---------- image ops ----------
C.prepImage = function (d, w, h) {
  const g = new Float32Array(w * h), b = new Float32Array(w * h), b2 = new Float32Array(w * h);
  for (let i = 0, j = 0; i < w * h; i++, j += 4) g[i] = 0.299 * d[j] + 0.587 * d[j + 1] + 0.114 * d[j + 2];
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) { const i = y * w + x, l = x > 0 ? i - 1 : i, r = x < w - 1 ? i + 1 : i; b[i] = 0.25 * g[l] + 0.5 * g[i] + 0.25 * g[r]; }
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) { const i = y * w + x, u = y > 0 ? i - w : i, dn = y < h - 1 ? i + w : i; b2[i] = 0.25 * b[u] + 0.5 * b[i] + 0.25 * b[dn]; }
  const gx = new Float32Array(w * h), gy = new Float32Array(w * h);
  for (let y = 1; y < h - 1; y++) for (let x = 1; x < w - 1; x++) { const i = y * w + x; gx[i] = 0.5 * (b2[i + 1] - b2[i - 1]); gy[i] = 0.5 * (b2[i + w] - b2[i - w]); }
  return { gray: b2, gx, gy };
};
const bil = (A, w, h, x, y) => { if (!(x >= 0 && y >= 0 && x <= w - 1 && y <= h - 1)) return 0; let x0 = Math.floor(x), y0 = Math.floor(y); if (x0 >= w - 1) x0 = w - 2; if (y0 >= h - 1) y0 = h - 2; const fx = x - x0, fy = y - y0, i = y0 * w + x0; return (A[i] * (1 - fx) + A[i + 1] * fx) * (1 - fy) + (A[i + w] * (1 - fx) + A[i + w + 1] * fx) * fy; };
C.cornerSubPix = function (img, p, win) {
  win = win || 5; if (!img.gx) return null; let q = p.slice(); const sg = 2 * (win * 0.6) ** 2;
  for (let it = 0; it < 40; it++) {
    let a = 0, b = 0, c = 0, bx = 0, by = 0;
    for (let dy = -win; dy <= win; dy++) for (let dx = -win; dx <= win; dx++) {
      const x = q[0] + dx, y = q[1] + dy, gx = bil(img.gx, img.w, img.h, x, y), gy = bil(img.gy, img.w, img.h, x, y), wt = Math.exp(-(dx * dx + dy * dy) / sg);
      const xx = wt * gx * gx, xy = wt * gx * gy, yy = wt * gy * gy; a += xx; b += xy; c += yy; bx += xx * x + xy * y; by += xy * x + yy * y;
    }
    const det = a * c - b * b, tr = a + c; if (!(tr > 0)) return null;
    const lmin = tr / 2 - Math.sqrt(tr * tr / 4 - det), lmax = tr / 2 + Math.sqrt(tr * tr / 4 - det);
    if (lmin / lmax < 0.08) return null;
    const nx = (c * bx - b * by) / det, ny = (a * by - b * bx) / det, sh = Math.hypot(nx - q[0], ny - q[1]); q = [nx, ny];
    if (Math.hypot(q[0] - p[0], q[1] - p[1]) > win * 1.2) return null; if (sh < 1e-3) break;
  }
  return q;
};
C.edgeAlong = function (img, p, n, range) {
  const st = 0.5, vals = []; let best = -1, bv = 0;
  for (let i = 0, s = -range; s <= range + 1e-9; s += st, i++) { const x = p[0] + s * n[0], y = p[1] + s * n[1], g = bil(img.gx, img.w, img.h, x, y) * n[0] + bil(img.gy, img.w, img.h, x, y) * n[1]; vals.push(g); if (Math.abs(g) > Math.abs(bv)) { bv = g; best = i; } }
  if (best <= 0 || best >= vals.length - 1) return null;
  const y0 = Math.abs(vals[best - 1]), y1 = Math.abs(vals[best]), y2 = Math.abs(vals[best + 1]), den = y0 - 2 * y1 + y2;
  const off = den < 0 ? Math.max(-0.5, Math.min(0.5, 0.5 * (y0 - y2) / den)) : 0, s = -range + (best + off) * st;
  return { p: [p[0] + s * n[0], p[1] + s * n[1]], g: bv };
};
// Fit a segment to an image edge. The fit is done in undistorted (ideal) space, so straight
// world edges stay straight under lens distortion; endpoints are returned in raw pixels.
C.refineSegment = function (img, a, b, cam) {
  if (!img.gx) return null;
  const id = cam ? (p => C.undistortPixel(cam, p[0], p[1])) : (p => p), rw = cam ? (p => C.distortPixel(cam, p[0], p[1])) : (p => p);
  let A = id(a), B = id(b); if (!isFinite(A[0]) || !isFinite(B[0])) return null;
  for (const range of [4, 1.5]) {
    const L = Math.hypot(B[0] - A[0], B[1] - A[1]); if (L < 6) return null;
    const K = Math.max(8, Math.min(240, Math.round(L / 3))), S = [];
    for (let i = 0; i < K; i++) {
      const t = 0.06 + 0.88 * i / (K - 1), q = [A[0] + t * (B[0] - A[0]), A[1] + t * (B[1] - A[1])], dt = 0.5 / L;
      const r0 = rw(q), r1 = rw([q[0] + dt * (B[0] - A[0]), q[1] + dt * (B[1] - A[1])]), tl = Math.hypot(r1[0] - r0[0], r1[1] - r0[1]); if (!(tl > 0)) continue;
      const n = [-(r1[1] - r0[1]) / tl, (r1[0] - r0[0]) / tl], ed = C.edgeAlong(img, r0, n, range);
      if (ed) { const u = id(ed.p); if (isFinite(u[0])) S.push({ p: u, g: ed.g }); }
    }
    if (S.length < 5) return null;
    const pos = S.filter(x => x.g > 0).length, sgn = pos >= S.length / 2 ? 1 : -1, mags = S.filter(x => x.g * sgn > 0).map(x => Math.abs(x.g)).sort((x, y) => x - y);
    const med = mags[Math.floor(mags.length / 2)] || 0; if (med < 1.5) return null;
    const P = S.filter(x => x.g * sgn > 0.4 * med).map(x => x.p); if (P.length < 5) return null;
    const ft = C.lineFitTLS(P), proj = q => { const t = (q[0] - ft.c[0]) * ft.dir[0] + (q[1] - ft.c[1]) * ft.dir[1]; return [ft.c[0] + t * ft.dir[0], ft.c[1] + t * ft.dir[1]]; };
    A = proj(A); B = proj(B);
  }
  const ra = rw(A), rb = rw(B);
  if (Math.hypot(ra[0] - a[0], ra[1] - a[1]) > 6 || Math.hypot(rb[0] - b[0], rb[1] - b[1]) > 6) return null;
  return [ra, rb];
};
C.refinePoly = function (img, pts) {
  if (!img.gx || pts.length < 2) return pts;
  return pts.map((p, i) => { const a = pts[Math.max(0, i - 1)], b = pts[Math.min(pts.length - 1, i + 1)], L = Math.hypot(b[0] - a[0], b[1] - a[1]); if (L < 1e-6) return p; const n = [-(b[1] - a[1]) / L, (b[0] - a[0]) / L], e = C.edgeAlong(img, p, n, 4); return e && Math.abs(e.g) > 1.5 ? e.p : p; });
};
C.undistortImage = function (d, w, h, cam) {
  const out = new Uint8ClampedArray(w * h * 4); let px = 0, py = 0;
  for (let v = 0; v < h; v++) for (let u = 0; u < w; u++) {
    const q = C.distortPixel(cam, u, v), x = q[0], y = q[1], o = (v * w + u) * 4;
    if (!(x >= 0 && y >= 0 && x <= w - 1 && y <= h - 1)) { out[o] = 11; out[o + 1] = 13; out[o + 2] = 15; out[o + 3] = 255; continue; }
    let x0 = Math.floor(x), y0 = Math.floor(y); if (x0 >= w - 1) x0 = w - 2; if (y0 >= h - 1) y0 = h - 2;
    const fx = x - x0, fy = y - y0, i = (y0 * w + x0) * 4, j = i + w * 4;
    for (let k = 0; k < 3; k++) out[o + k] = (d[i + k] * (1 - fx) + d[i + 4 + k] * fx) * (1 - fy) + (d[j + k] * (1 - fx) + d[j + 4 + k] * fx) * fy;
    out[o + 3] = 255;
  }
  return out;
};
// ray-traced room corner: floor Z=0 (0..8)², walls X=0 and Y=0 (height 5), 1-unit checkers
C.renderSynthetic = function (cam) {
  const w = cam.w, h = cam.h, out = new Uint8ClampedArray(w * h * 4), R = cam.R, O = C.camCenter(cam), rnd = rng(7);
  const r00 = R[0][0], r01 = R[0][1], r02 = R[0][2], r10 = R[1][0], r11 = R[1][1], r12 = R[1][2], r20 = R[2][0], r21 = R[2][1], r22 = R[2][2];
  const O0 = O[0], O1 = O[1], O2 = O[2], k1 = cam.k1, k2 = cam.k2, k3 = cam.k3, p1 = cam.p1, p2 = cam.p2;
  const PAL = [[38, 40, 44], [232, 230, 224], [52, 64, 82], [206, 214, 224], [86, 66, 52], [224, 214, 200], [26, 29, 33]];
  const SX = [-0.25, 0.25, -0.25, 0.25], SY = [-0.25, -0.25, 0.25, 0.25];
  for (let v = 0; v < h; v++) {
    let x = (v - cam.cy) / cam.fy * 0, y = 0, warm = false;
    for (let u = 0; u < w; u++) {
      let ar = 0, ag = 0, ab = 0;
      for (let s = 0; s < 4; s++) {
        const xd = (u + SX[s] - cam.cx) / cam.fx, yd = (v + SY[s] - cam.cy) / cam.fy;
        if (!warm) { x = xd; y = yd; warm = true; }
        for (let i = 0; i < 12; i++) {
          const r2 = x * x + y * y, rad = 1 + r2 * (k1 + r2 * (k2 + r2 * k3)), Rr = k1 + r2 * (2 * k2 + 3 * k3 * r2);
          const ex = x * rad + 2 * p1 * x * y + p2 * (r2 + 2 * x * x) - xd, ey = y * rad + p1 * (r2 + 2 * y * y) + 2 * p2 * x * y - yd;
          const a = rad + 2 * x * x * Rr + 2 * p1 * y + 6 * p2 * x, b = 2 * x * y * Rr + 2 * p1 * x + 2 * p2 * y, d = rad + 2 * y * y * Rr + 6 * p1 * y + 2 * p2 * x, det = a * d - b * b;
          const dx = (d * ex - b * ey) / det, dy = (a * ey - b * ex) / det; x -= dx; y -= dy; if (dx * dx + dy * dy < 1e-26) break;
        }
        // world ray = R^T [x,y,1]
        const d0 = r00 * x + r10 * y + r20, d1 = r01 * x + r11 * y + r21, d2 = r02 * x + r12 * y + r22;
        let best = Infinity, col = 6, t, P0, P1, P2;
        if (Math.abs(d2) > 1e-12) { t = -O2 / d2; if (t > 0 && t < best) { P0 = O0 + t * d0; P1 = O1 + t * d1; if (P0 >= 0 && P0 <= 8 && P1 >= 0 && P1 <= 8) { best = t; col = ((Math.floor(P0) + Math.floor(P1)) & 1) ? 1 : 0; } } }
        if (Math.abs(d0) > 1e-12) { t = -O0 / d0; if (t > 0 && t < best) { P1 = O1 + t * d1; P2 = O2 + t * d2; if (P1 >= 0 && P1 <= 8 && P2 >= 0 && P2 <= 5) { best = t; col = ((Math.floor(P1) + Math.floor(P2)) & 1) ? 3 : 2; } } }
        if (Math.abs(d1) > 1e-12) { t = -O1 / d1; if (t > 0 && t < best) { P0 = O0 + t * d0; P2 = O2 + t * d2; if (P0 >= 0 && P0 <= 8 && P2 >= 0 && P2 <= 5) { best = t; col = ((Math.floor(P0) + Math.floor(P2)) & 1) ? 5 : 4; } } }
        const c = PAL[col]; ar += c[0]; ag += c[1]; ab += c[2];
      }
      const o = (v * w + u) * 4, nz = (rnd() - 0.5) * 4; out[o] = ar / 4 + nz; out[o + 1] = ag / 4 + nz; out[o + 2] = ab / 4 + nz; out[o + 3] = 255;
    }
  }
  return out;
};
if (typeof window !== 'undefined') window.Calib = C;
if (typeof module !== 'undefined') module.exports = C;
})();
