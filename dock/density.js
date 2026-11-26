/* Electron density: read a CCP4/MRC map and compute the real-space R value (RSR) of a set of atoms against it.
   RSR = sum|rho_obs - s*rho_calc| / sum|rho_obs + s*rho_calc| over the grid points near the ligand, with s fitted by least squares.
   rho_calc is a sum of Gaussians (weight = atomic number, width from B-factor and resolution). Lower is better; the crystal ligand itself is the reference. */
(function (root) {
'use strict';
const Z = { H: 1, C: 6, N: 7, O: 8, F: 9, P: 15, S: 16, Cl: 17, Br: 35, I: 53, Se: 34, Zn: 30, Fe: 26, Mg: 12, Ca: 20, Na: 11, K: 19, Mn: 25, Cu: 29, Co: 27, Ni: 28, B: 5, Si: 14 };
const RAD = Math.PI / 180;

function cellMatrix(c) {
  const [a, b, cc] = c, ca = Math.cos(c[3] * RAD), cb = Math.cos(c[4] * RAD), cg = Math.cos(c[5] * RAD), sg = Math.sin(c[5] * RAD);
  const v = Math.sqrt(Math.max(1e-9, 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg));
  /* fractional -> Cartesian (PDB convention), upper triangular */
  return { m00: a, m01: b * cg, m02: cc * cb, m11: b * sg, m12: cc * (ca - cb * cg) / sg, m22: cc * v / sg };
}
function toFrac(M, p) {
  const z = p[2] / M.m22, y = (p[1] - M.m12 * z) / M.m11, x = (p[0] - M.m01 * y - M.m02 * z) / M.m00; return [x, y, z];
}
function toCart(M, f) { return [M.m00 * f[0] + M.m01 * f[1] + M.m02 * f[2], M.m11 * f[1] + M.m12 * f[2], M.m22 * f[2]]; }

/* CCP4 / MRC2014 map -> grid object */
function parseMap(buf) {
  if (buf.byteLength < 1024) throw new Error('That file is too small to be a CCP4 or MRC map.');
  const dv = new DataView(buf); let le = true;
  const ms = dv.getUint8(212); if (ms === 0x11) le = false;
  let mode = dv.getInt32(12, le); if (!(mode >= 0 && mode <= 6)) { le = !le; mode = dv.getInt32(12, le); }
  if (!(mode >= 0 && mode <= 6)) throw new Error('Not a CCP4 or MRC map (unreadable header).');
  const I = o => dv.getInt32(o, le), F = o => dv.getFloat32(o, le);
  const dims = [I(0), I(4), I(8)], start = [I(16), I(20), I(24)], m = [I(28), I(32), I(36)];
  const cell = [F(40), F(44), F(48), F(52), F(56), F(60)], axis = [I(64) - 1, I(68) - 1, I(72) - 1], nsym = I(92);
  if (dims.some(v => v <= 0) || m.some(v => v <= 0) || cell.slice(0, 3).some(v => !(v > 0))) throw new Error('The map header is not valid.');
  if ([...axis].sort().join() !== '0,1,2') throw new Error('Unsupported map axis order.');
  const n = dims[0] * dims[1] * dims[2], off = 1024 + nsym, data = new Float32Array(n);
  const sz = { 0: 1, 1: 2, 2: 4, 6: 2 }[mode]; if (sz === undefined) throw new Error('Unsupported map data type (mode ' + mode + ').');
  if (buf.byteLength < off + n * sz) throw new Error('The map file is truncated.');
  if (mode === 2) { if (le && (off % 4 === 0)) data.set(new Float32Array(buf, off, n)); else for (let i = 0; i < n; i++) data[i] = dv.getFloat32(off + 4 * i, le); }
  else for (let i = 0; i < n; i++) data[i] = mode === 0 ? dv.getInt8(off + i) : mode === 1 ? dv.getInt16(off + 2 * i, le) : dv.getUint16(off + 2 * i, le);
  let sum = 0, sq = 0; for (let i = 0; i < n; i++) { sum += data[i]; sq += data[i] * data[i]; }
  const mean = sum / n, sd = Math.sqrt(Math.max(0, sq / n - mean * mean));
  /* which file dimension carries crystal axis k */
  const dim = [0, 0, 0]; for (let d = 0; d < 3; d++) dim[axis[d]] = d;
  return { dims, start, m, cell, axis, dim, data, M: cellMatrix(cell), mean, sd, n };
}
/* value at integer crystal grid index (i,j,k), wrapped by the cell; null when the map does not cover it */
function nodeValue(g, ijk) {
  const p = [0, 0, 0];
  for (let k = 0; k < 3; k++) { const d = g.dim[k]; let v = ((ijk[k] % g.m[k]) + g.m[k]) % g.m[k] - g.start[d]; if (v < 0) v += g.m[k]; if (v < 0 || v >= g.dims[d]) return null; p[d] = v; }
  return g.data[p[0] + g.dims[0] * (p[1] + g.dims[1] * p[2])];
}
/* grid nodes (unwrapped integer indices) within r of any of the atoms: returns Map key -> {ijk, xyz} */
function nodesNear(g, atoms, r, into) {
  const out = into || new Map(), M = g.M;
  /* fractional half-widths of a sphere of radius r: r*|row of inverse| ; use generous bound via cell lengths */
  const hw = [r / (g.cell[0] * 0.5) * 0.5, r / (g.cell[1] * 0.5) * 0.5, r / (g.cell[2] * 0.5) * 0.5];
  const sa = Math.sin(g.cell[3] * RAD), sb = Math.sin(g.cell[4] * RAD), sg = Math.sin(g.cell[5] * RAD);
  const w = [r / (g.cell[0] * sb * sg), r / (g.cell[1] * sa * sg), r / (g.cell[2] * sa * sb)].map(v => Math.min(v * 1.05 + 1e-6, 1));
  void hw; const r2 = r * r;
  for (const a of atoms) {
    const f = toFrac(M, a.xyz);
    const lo = [0, 1, 2].map(k => Math.floor((f[k] - w[k]) * g.m[k])), hi = [0, 1, 2].map(k => Math.ceil((f[k] + w[k]) * g.m[k]));
    for (let i = lo[0]; i <= hi[0]; i++) for (let j = lo[1]; j <= hi[1]; j++) for (let k = lo[2]; k <= hi[2]; k++) {
      const key = i + ',' + j + ',' + k; if (out.has(key)) continue;
      const c = toCart(M, [i / g.m[0], j / g.m[1], k / g.m[2]]), dx = c[0] - a.xyz[0], dy = c[1] - a.xyz[1], dz = c[2] - a.xyz[2];
      if (dx * dx + dy * dy + dz * dz <= r2) out.set(key, { ijk: [i, j, k], xyz: c });
    }
  }
  return out;
}
function sigmaOf(B, res) { return Math.sqrt(Math.max(5, B) / (8 * Math.PI * Math.PI) + Math.pow(0.25 * (res || 2), 2)); }

/* atoms: [{xyz:[x,y,z], el}] for the pose; ref: the crystal ligand atoms (defines the common region); opt: {res, B, r} */
function rsr(g, atoms, ref, opt = {}) {
  const sig = sigmaOf(opt.B || 20, opt.res || 2), rr = opt.r || 2.0, cut = 4 * sig;
  const region = nodesNear(g, ref || atoms, rr + 0.5); if (ref) nodesNear(g, atoms, rr, region);
  const obs = [], calc = [];
  for (const nd of region.values()) {
    const v = nodeValue(g, nd.ijk); if (v === null) continue;
    let c = 0; for (const a of atoms) { const dx = nd.xyz[0] - a.xyz[0], dy = nd.xyz[1] - a.xyz[1], dz = nd.xyz[2] - a.xyz[2], d2 = dx * dx + dy * dy + dz * dz; if (d2 < cut * cut) c += (Z[a.el] || 6) * Math.exp(-d2 / (2 * sig * sig)); }
    obs.push(v); calc.push(c);
  }
  if (obs.length < 20) return { rsr: NaN, n: obs.length, scale: 0 };
  let so = 0, sc = 0; for (let i = 0; i < obs.length; i++) { so += obs[i] * calc[i]; sc += calc[i] * calc[i]; }
  const s = sc > 0 ? so / sc : 0; let num = 0, den = 0;
  for (let i = 0; i < obs.length; i++) { num += Math.abs(obs[i] - s * calc[i]); den += Math.abs(obs[i] + s * calc[i]); }
  return { rsr: den > 0 ? num / den : NaN, n: obs.length, scale: s };
}
/* tiny synthetic map writer (for tests): Gaussian density of atoms on an orthogonal cell */
function synthMap(atoms, cell, step, sig, noise, rnd) {
  const m = cell.slice(0, 3).map(v => Math.round(v / step)), n = m[0] * m[1] * m[2], buf = new ArrayBuffer(1024 + 4 * n), dv = new DataView(buf);
  const w = (o, v) => dv.setInt32(o, v, true); w(0, m[0]); w(4, m[1]); w(8, m[2]); w(12, 2); w(16, 0); w(20, 0); w(24, 0); w(28, m[0]); w(32, m[1]); w(36, m[2]);
  cell.forEach((v, i) => dv.setFloat32(40 + 4 * i, v, true)); w(64, 1); w(68, 2); w(72, 3); dv.setUint8(208, 0x4d); dv.setUint8(209, 0x41); dv.setUint8(210, 0x50); dv.setUint8(212, 0x44); dv.setUint8(213, 0x44);
  const f = new Float32Array(buf, 1024, n), cut = 4 * sig;
  for (const a of atoms) { const lo = [0, 1, 2].map(k => Math.floor((a.xyz[k] - cut) / step)), hi = [0, 1, 2].map(k => Math.ceil((a.xyz[k] + cut) / step));
    for (let i = lo[0]; i <= hi[0]; i++) for (let j = lo[1]; j <= hi[1]; j++) for (let k = lo[2]; k <= hi[2]; k++) { const dx = i * step - a.xyz[0], dy = j * step - a.xyz[1], dz = k * step - a.xyz[2], d2 = dx * dx + dy * dy + dz * dz; if (d2 > cut * cut) continue;
      const I = ((i % m[0]) + m[0]) % m[0], J = ((j % m[1]) + m[1]) % m[1], K = ((k % m[2]) + m[2]) % m[2]; f[I + m[0] * (J + m[1] * K)] += (Z[a.el] || 6) * Math.exp(-d2 / (2 * sig * sig)); } }
  if (noise) for (let i = 0; i < n; i++) f[i] += noise * (rnd() - 0.5);
  return buf;
}
const API = { parseMap, rsr, nodeValue, nodesNear, sigmaOf, synthMap, cellMatrix, toFrac, toCart };
if (typeof module !== 'undefined' && module.exports) module.exports = API; else root.Density = API;
})(typeof window !== 'undefined' ? window : globalThis);
