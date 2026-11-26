/* ============================ nmr.js ============================
   CompChemX NMR: an empirical NMR predictor and spectrum simulator. No libraries, no network.

   Input : a V2000 molfile (from CompChemX Draw).
   Output: predicted 1H and 13C shifts, first-order multiplets and a simulated spectrum.

   Method (all additive increment schemes, the kind found in Pretsch, "Structure Determination of Organic Compounds"):
     1H  alkyl C-H   base value + alpha and beta substituent increments (Shoolery / Curphy-Morrison style)
         alkene C-H  5.25 + gem, cis and trans increments (Pascual-Meier-Simon)
         aromatic    7.27 + ortho, meta and para increments; base values for pyridine-type and five-ring heteroaromatics
         aldehyde, alkyne, O-H, N-H, S-H  typical values for the class
     13C alkyl       -2.3 + Grant-Paul counts, branching corrections and substituent alpha / beta / gamma increments
         alkene      123.3 + increments;  aromatic 128.5 + increments;  C=O, C#N and C#C by class
   Equivalent atoms (same topological environment) share one signal. Multiplets are first order (n+1 rule with
   typical 3J values: 7 Hz aliphatic, 7.8 Hz aromatic ortho, 1.8 Hz meta, cis 10 / trans 17 Hz alkene).

   What it is not: a quantum or database predictor. Expect about +/-0.3 ppm for 1H and +/-5 ppm for 13C on ordinary
   organic molecules. It does not model diastereotopic splitting, second-order effects, conformers or solvent shifts of
   C-H protons.                                                                                                    */
(function (root) {
'use strict';

/* ------------------------------------------------------------------ molfile and graph */
function parseMolfile(text) {
  const ln = String(text).replace(/\r/g, '').split('\n');
  if (!/V2000/.test(ln[3] || '')) throw new Error('Expected a V2000 molfile.');
  const na = parseInt(ln[3].substr(0, 3)), nb = parseInt(ln[3].substr(3, 3));
  if (!(na > 0)) throw new Error('The drawing has no atoms.');
  const atoms = [], bonds = [];
  const CHGCODE = { 1: 3, 2: 2, 3: 1, 5: -1, 6: -2, 7: -3 };
  for (let i = 0; i < na; i++) {
    const l = ln[4 + i];
    atoms.push({ el: l.substr(31, 3).trim(), x: parseFloat(l.substr(0, 10)) || 0, y: parseFloat(l.substr(10, 10)) || 0, chg: CHGCODE[parseInt(l.substr(36, 3))] || 0 });
  }
  for (let i = 0; i < nb; i++) {
    const l = ln[4 + na + i];
    bonds.push({ a: parseInt(l.substr(0, 3)) - 1, b: parseInt(l.substr(3, 3)) - 1, o: parseInt(l.substr(6, 3)) || 1 });
  }
  let sawM = false;
  for (const l of ln) if (/^M  CHG/.test(l)) {
    if (!sawM) { atoms.forEach(a => a.chg = 0); sawM = true; }
    const n = parseInt(l.substr(6, 3)); for (let i = 0; i < n; i++) { const k = parseInt(l.substr(10 + i * 8, 4)) - 1, v = parseInt(l.substr(14 + i * 8, 4)); if (atoms[k]) atoms[k].chg = v; }
  }
  return { atoms, bonds };
}

const ALLOWED = { C: [4], N: [3], O: [2], S: [2, 4, 6], P: [3, 5], F: [1], Cl: [1], Br: [1], I: [1], B: [3], Si: [4], Se: [2, 4, 6] };
function allowedValences(el, chg) {
  const base = ALLOWED[el] || [4];
  if (!chg) return base;
  if (el === 'N' || el === 'P') return [base[0] + chg];
  if (el === 'O' || el === 'S') return [base[0] + chg];
  if (el === 'C') return [3];
  if (el === 'B') return [chg < 0 ? 4 : 2];
  return [Math.max(0, base[0] + chg)];
}

/* Build the heavy-atom graph: explicit hydrogens are folded into counts, hydrogen counts come from valence. */
function build(mol) {
  const idx = new Array(mol.atoms.length).fill(-1), atoms = [];
  mol.atoms.forEach((a, i) => { if (!/^[HDT]$/.test(a.el)) { idx[i] = atoms.length; atoms.push({ el: a.el, chg: a.chg || 0, x: a.x, y: a.y, nb: [], src: i, nH: 0, arom: false, rings: [] }); } });
  const bonds = [];
  for (const b of mol.bonds) {
    const i = idx[b.a], j = idx[b.b]; if (i < 0 || j < 0 || i === j) continue;
    const bd = { a: i, b: j, o: Math.min(3, Math.max(1, b.o)), ar: false, i: bonds.length };
    bonds.push(bd); atoms[i].nb.push({ j, bd }); atoms[j].nb.push({ j: i, bd });
  }
  const G = { atoms, bonds, n: atoms.length, unknown: [] };
  for (const a of atoms) {
    if (!ALLOWED[a.el]) G.unknown.push(a.el);
    let bo = 0; for (const { bd } of a.nb) bo += bd.o;
    const v = allowedValences(a.el, a.chg).find(x => x >= bo);
    a.nH = v === undefined ? 0 : v - bo;
  }
  findRings(G); perceiveAromatic(G); symmetry(G);
  return G;
}

/* ------------------------------------------------------------------ rings and aromaticity */
function findRings(G, maxSize = 8) {
  const n = G.n, seen = new Set(), cycles = [];
  for (let s = 0; s < n && cycles.length < 4000; s++) {
    const path = [s], onPath = new Array(n).fill(false); onPath[s] = true;
    (function dfs(u) {
      if (cycles.length >= 4000) return;
      for (const { j } of G.atoms[u].nb) {
        if (j === s && path.length >= 3) { const key = path.slice().sort((a, b) => a - b).join(','); if (!seen.has(key)) { seen.add(key); cycles.push(path.slice()); } continue; }
        if (j <= s || onPath[j] || path.length >= maxSize) continue;
        path.push(j); onPath[j] = true; dfs(j); path.pop(); onPath[j] = false;
      }
    })(s);
  }
  cycles.sort((a, b) => a.length - b.length);
  const bidx = new Map(); G.bonds.forEach(b => { bidx.set(b.a * n + b.b, b.i); bidx.set(b.b * n + b.a, b.i); });
  const comps = (() => { const c = new Array(n).fill(-1); let k = 0; for (let i = 0; i < n; i++) if (c[i] < 0) { const q = [i]; c[i] = k; while (q.length) { const u = q.pop(); for (const { j } of G.atoms[u].nb) if (c[j] < 0) { c[j] = k; q.push(j); } } k++; } return k; })();
  const want = G.bonds.length - n + comps, basis = [], rings = [];
  for (const cyc of cycles) {
    if (rings.length >= want) break;
    const v = new Uint8Array(G.bonds.length);
    for (let i = 0; i < cyc.length; i++) v[bidx.get(cyc[i] * n + cyc[(i + 1) % cyc.length])] ^= 1;
    const w = v.slice();
    for (const b of basis) if (w[b.p]) for (let k = 0; k < w.length; k++) w[k] ^= b.v[k];
    const p = w.indexOf(1); if (p < 0) continue;
    basis.push({ p, v: w }); rings.push({ atoms: cyc, size: cyc.length, bonds: cyc.map((x, i) => bidx.get(x * n + cyc[(i + 1) % cyc.length])), arom: false });
  }
  G.rings = rings;
  const inRing = new Array(n).fill(false); rings.forEach(r => r.atoms.forEach(a => inRing[a] = true));
  G.inRing = inRing;
}

function perceiveAromatic(G) {
  for (const r of G.rings) {
    if (r.size < 5 || r.size > 7) continue;
    const set = new Set(r.atoms); let e = 0, ok = true;
    for (const i of r.atoms) {
      const a = G.atoms[i];
      let dbl = false, trip = false;
      for (const { j, bd } of a.nb) { if (bd.o === 2 && (G.inRing[j]) && !(a.el === 'S' && G.atoms[j].el === 'O')) dbl = true; if (bd.o === 3) trip = true; }
      if (trip) { ok = false; break; }
      if (dbl) e += 1;
      else if (a.nb.some(x => x.bd.o === 2)) e += 0;                                  /* exocyclic C=O etc: contributes nothing */
      else if (a.el === 'O' || a.el === 'S' || a.el === 'Se') e += 2;
      else if (a.el === 'N' && a.chg <= 0) e += 2;
      else if (a.el === 'C' && a.chg < 0) e += 2;
      else if (a.el === 'C' && a.chg > 0) e += 0;
      else if (a.el === 'B') e += 0;
      else { ok = false; break; }
    }
    if (ok && (e === 6 || e === 2 || e === 10)) r.arom = true;
  }
  G.rings.forEach((r, k) => { if (r.arom) { r.atoms.forEach(i => { G.atoms[i].arom = true; G.atoms[i].rings.push(k); }); r.bonds.forEach(b => G.bonds[b].ar = true); } });
  G.atoms.forEach((a, i) => { G.rings.forEach((r, k) => { if (!r.arom && r.atoms.includes(i)) a.nrings = (a.nrings || []).concat(k); }); });
}

/* ------------------------------------------------------------------ symmetry classes (Morgan style refinement) */
function symmetry(G) {
  const n = G.n; let cls = new Array(n);
  const init = G.atoms.map(a => [a.el, a.chg, a.nH, a.arom ? 1 : 0, a.nb.length].join('/'));
  let map = {}; init.forEach(k => { if (!(k in map)) map[k] = Object.keys(map).length; }); cls = init.map(k => map[k]);
  for (let it = 0; it < n + 2; it++) {
    const keys = G.atoms.map((a, i) => cls[i] + '|' + a.nb.map(({ j, bd }) => (bd.ar ? 'a' : bd.o) + ':' + cls[j]).sort().join(','));
    map = {}; const nc = keys.map(k => { if (!(k in map)) map[k] = Object.keys(map).length; return map[k]; });
    const before = new Set(cls).size, after = Object.keys(map).length; cls = nc; if (after === before) break;
  }
  G.atoms.forEach((a, i) => a.cls = cls[i]);
}

/* ------------------------------------------------------------------ chemistry helpers */
const at = (G, i) => G.atoms[i];
const heavyN = (G, i) => G.atoms[i].nb;
const hasPi = (G, i) => G.atoms[i].nb.some(x => x.bd.o > 1) || G.atoms[i].arom;
const isSp3C = (G, i) => { const a = G.atoms[i]; return a.el === 'C' && !a.arom && !a.nb.some(x => x.bd.o > 1) && a.chg === 0; };
const dblO = (G, i) => { const a = G.atoms[i]; for (const { j, bd } of a.nb) if (bd.o === 2 && G.atoms[j].el === 'O') return j; return -1; };
const isCarbonyl = (G, i) => G.atoms[i].el === 'C' && !G.atoms[i].arom && dblO(G, i) >= 0;
const isAlkeneC = (G, i) => { const a = G.atoms[i]; return a.el === 'C' && !a.arom && a.nb.some(x => x.bd.o === 2 && G.atoms[x.j].el === 'C' && !G.atoms[x.j].arom); };
const tripleNbr = (G, i) => { for (const { j, bd } of G.atoms[i].nb) if (bd.o === 3) return j; return -1; };
const isAcylN = (G, n, except) => G.atoms[n].nb.some(x => x.j !== except && (isCarbonyl(G, x.j) || (G.atoms[x.j].el === 'S' && G.atoms[x.j].nb.some(y => y.bd.o === 2 && G.atoms[y.j].el === 'O'))));
const hasArNbr = (G, n, except) => G.atoms[n].nb.some(x => x.j !== except && G.atoms[x.j].arom);

/* Which kind of group is atom s, seen from atom c that it is attached to? */
function subKey(G, c, s) {
  const a = G.atoms[s], el = a.el, others = a.nb.filter(x => x.j !== c);
  if (el === 'F' || el === 'Cl' || el === 'Br' || el === 'I') return el;
  if (el === 'C') {
    if (a.arom) return 'Ph';
    const t = tripleNbr(G, s);
    if (t >= 0) return G.atoms[t].el === 'N' ? 'CN' : 'alkynyl';
    const o = dblO(G, s);
    if (o >= 0) {
      const rest = others.filter(x => x.j !== o);
      if (a.nH > 0 && rest.length === 0) return 'CHO';
      const oo = rest.find(x => G.atoms[x.j].el === 'O'); if (oo) return G.atoms[oo.j].nH > 0 || G.atoms[oo.j].chg < 0 ? 'COOH' : 'COOR';
      if (rest.some(x => G.atoms[x.j].el === 'N')) return 'CON';
      if (rest.some(x => G.atoms[x.j].el === 'Cl')) return 'COCl';
      return 'COR';
    }
    if (a.nb.some(x => x.bd.o === 2)) return 'vinyl';
    const hs = others.map(x => G.atoms[x.j].el);
    if (hs.filter(e => e === 'F').length >= 3) return 'CF3';
    if (hs.includes('O')) return 'CH2O';
    if (hs.includes('N')) return 'CH2N';
    if (hs.some(e => e === 'Cl' || e === 'Br' || e === 'I')) return 'CH2X';
    const nC = others.length;
    return nC === 0 ? 'CH3' : nC === 1 ? 'Et' : nC === 2 ? 'iPr' : 'tBu';
  }
  if (el === 'N') {
    if (a.arom) return 'NArom';
    if (others.filter(x => G.atoms[x.j].el === 'O').length >= 2 || (a.chg > 0 && others.some(x => G.atoms[x.j].el === 'O'))) return 'NO2';
    if (a.chg > 0 || a.nb.length + a.nH >= 4) return 'NR3';
    if (isAcylN(G, s, c)) return 'NHAc';
    return a.nH === 2 ? 'NH2' : a.nH === 1 ? 'NHR' : 'NR2';
  }
  if (el === 'O') {
    if (others.some(x => isCarbonyl(G, x.j) || (G.atoms[x.j].el === 'S' && G.atoms[x.j].nb.some(y => y.bd.o === 2)) || G.atoms[x.j].el === 'P')) return 'OAc';
    if (others.some(x => G.atoms[x.j].arom)) return 'OPh';
    return a.nH > 0 ? 'OH' : 'OR';
  }
  if (el === 'S') {
    const no = a.nb.filter(x => x.bd.o === 2 && G.atoms[x.j].el === 'O').length;
    return no >= 2 ? 'SO2' : no === 1 ? 'SO' : a.nH > 0 ? 'SH' : 'SR';
  }
  if (el === 'P') return 'P'; if (el === 'Si') return 'Si'; if (el === 'B') return 'B';
  return 'X';
}

const RING_BASE_CH2 = { 3: 0.25, 4: 1.95, 5: 1.50, 6: 1.43, 7: 1.50 }, RING_BASE_CH = { 3: 0.80, 4: 1.95, 5: 1.65, 6: 1.55, 7: 1.60 };
function smallestRing(G, i) { let s = 99; G.rings.forEach(r => { if (!r.arom && r.atoms.includes(i) && r.size < s) s = r.size; }); return s === 99 ? 0 : s; }

/* ------------------------------------------------------------------ 1H increment tables */
const H_ALPHA = { F: 3.2, Cl: 2.3, Br: 2.1, I: 1.7, OH: 2.35, OR: 2.35, OPh: 2.9, OAc: 2.95, NH2: 1.45, NHR: 1.45, NR2: 1.45, NHAc: 2.0, NR3: 2.3, NO2: 3.45, NArom: 2.7, SH: 1.25, SR: 1.25, SO: 1.7, SO2: 2.0,
  Ph: 1.45, vinyl: 0.8, alkynyl: 0.9, CN: 1.15, CHO: 1.35, COR: 1.3, COOH: 1.15, COOR: 1.15, CON: 1.1, COCl: 1.5, P: 0.9, Si: -0.7, B: 0, X: 0 };
const H_BETA = { F: 0.3, Cl: 0.6, Br: 0.8, I: 0.95, OH: 0.35, OR: 0.35, OPh: 0.4, OAc: 0.45, NH2: 0.25, NHR: 0.25, NR2: 0.25, NHAc: 0.3, NR3: 0.4, NO2: 0.65, NArom: 0.4, SH: 0.3, SR: 0.3, SO: 0.4, SO2: 0.4,
  Ph: 0.35, vinyl: 0.2, alkynyl: 0.2, CN: 0.25, CHO: 0.2, COR: 0.2, COOH: 0.2, COOR: 0.2, CON: 0.2, COCl: 0.25, CF3: 0.7, P: 0.2 };
/* alkene: [gem, cis, trans] */
const ALK_ALKYL = [0.45, -0.22, -0.28];
const H_VINYL = { CH3: ALK_ALKYL, Et: ALK_ALKYL, iPr: ALK_ALKYL, tBu: [0.55, -0.21, -0.28], CH2O: [0.64, -0.01, -0.02], CH2X: [0.70, 0.11, -0.04], CH2N: [0.58, -0.10, -0.08], CF3: [0.66, 0.61, 0.32],
  vinyl: [1.00, -0.09, -0.23], alkynyl: [0.47, 0.38, 0.12], Ph: [1.38, 0.36, -0.07], CHO: [1.03, 0.97, 1.21], COR: [1.10, 1.12, 0.87], COOH: [1.00, 1.35, 0.74], COOR: [0.84, 1.15, 0.56], CON: [1.37, 0.98, 0.46],
  COCl: [1.10, 1.46, 1.01], CN: [0.27, 0.75, 0.55], OH: [1.22, -1.07, -1.21], OR: [1.22, -1.07, -1.21], OPh: [1.21, -0.60, -1.00], OAc: [2.11, -0.35, -0.64], NH2: [0.80, -1.26, -1.21], NHR: [0.80, -1.26, -1.21],
  NR2: [0.80, -1.26, -1.21], NHAc: [2.08, -0.57, -0.72], NArom: [2.08, -0.57, -0.72], NO2: [1.87, 1.30, 0.62], F: [1.54, -0.40, -1.02], Cl: [1.08, 0.18, 0.13], Br: [1.07, 0.45, 0.55], I: [1.14, 0.81, 0.88],
  SH: [1.11, -0.29, -0.13], SR: [1.11, -0.29, -0.13], SO2: [1.55, 1.16, 0.93], SO: [1.55, 1.16, 0.93], Si: [0.0, 0.0, 0.0] };
/* aromatic: [ortho, meta, para] for 1H and [ipso, ortho, meta, para] for 13C */
const AR_H = { CH3: [-0.20, -0.12, -0.22], Et: [-0.14, -0.06, -0.17], iPr: [-0.13, -0.08, -0.18], tBu: [0.02, -0.08, -0.21], CH2O: [-0.07, -0.07, -0.07], CH2X: [0.0, 0.0, 0.0], CH2N: [-0.05, -0.05, -0.08],
  CF3: [0.32, 0.14, 0.20], vinyl: [0.04, -0.04, -0.12], alkynyl: [0.15, -0.02, -0.01], Ph: [0.37, 0.20, 0.10], CHO: [0.56, 0.22, 0.29], COR: [0.62, 0.14, 0.21], COOH: [0.85, 0.18, 0.27], COOR: [0.71, 0.11, 0.21],
  CON: [0.61, 0.10, 0.17], COCl: [0.84, 0.22, 0.36], CN: [0.36, 0.18, 0.28], NO2: [0.95, 0.26, 0.38], OH: [-0.56, -0.12, -0.45], OR: [-0.48, -0.09, -0.44], OPh: [-0.29, -0.05, -0.23], OAc: [-0.25, 0.03, -0.13],
  NH2: [-0.75, -0.25, -0.65], NHR: [-0.80, -0.22, -0.68], NR2: [-0.66, -0.18, -0.67], NHAc: [0.12, -0.07, -0.28], NArom: [0.12, -0.07, -0.28], NR3: [0.69, 0.36, 0.31], F: [-0.26, 0.00, -0.20], Cl: [0.03, -0.02, -0.09],
  Br: [0.18, -0.08, -0.04], I: [0.39, -0.21, 0.00], SH: [-0.08, -0.16, -0.22], SR: [-0.08, -0.10, -0.24], SO: [0.45, 0.2, 0.2], SO2: [0.60, 0.26, 0.33], Si: [0.22, -0.02, -0.02], P: [0.1, 0.0, 0.0] };
const AR_C = { CH3: [9.3, 0.7, -0.1, -3.0], Et: [15.6, -0.5, 0.0, -2.6], iPr: [20.1, -2.0, 0.0, -2.5], tBu: [22.2, -3.4, -0.4, -3.1], CH2O: [12.4, -1.2, 0.2, -1.1], CH2X: [9.3, 0.2, 0.2, -0.2], CH2N: [14.0, -1.5, 0.0, -1.5],
  CF3: [2.6, -3.1, 0.4, 3.4], vinyl: [8.9, -2.3, -0.1, -0.8], alkynyl: [-6.2, 3.6, -0.4, -0.3], Ph: [13.0, -1.1, 0.5, -1.0], CHO: [8.2, 1.2, 0.5, 5.8], COR: [9.1, 0.1, 0.0, 4.2], COOH: [2.1, 1.6, -0.1, 5.2],
  COOR: [2.0, 1.2, -0.1, 4.3], CON: [5.0, -1.2, 0.0, 3.4], COCl: [4.6, 2.4, 0.3, 6.2], CN: [-15.7, 3.6, 0.7, 4.3], NO2: [19.9, -4.9, 0.9, 6.1], OH: [26.9, -12.7, 1.4, -7.3], OR: [31.4, -14.4, 1.0, -7.7],
  OPh: [27.6, -11.2, -0.3, -6.9], OAc: [22.4, -7.1, -0.4, -3.2], NH2: [18.2, -13.4, 0.8, -10.0], NHR: [21.7, -16.2, 0.7, -11.8], NR2: [22.4, -15.7, 0.8, -11.8], NHAc: [11.1, -9.9, 0.2, -5.6], NArom: [11.1, -9.9, 0.2, -5.6],
  NR3: [19.8, -5.3, 1.7, 2.4], F: [35.1, -14.3, 0.9, -4.4], Cl: [6.4, 0.2, 1.0, -1.9], Br: [-5.4, 3.3, 2.2, -1.0], I: [-32.3, 9.9, 2.6, -0.4], SH: [2.2, 0.7, 0.4, -3.3], SR: [10.0, -1.9, 0.2, -3.6],
  SO: [18.0, -3.0, 1.0, 2.0], SO2: [12.3, -1.4, 0.8, 5.1], Si: [11.0, 4.4, -0.8, 0.0], P: [9.0, 3.0, 0.0, 0.0] };
/* sp3 carbon substituent increments [alpha, beta, gamma] (Pretsch) */
const C_SUB = { F: [68, 9, -4], Cl: [31, 11, -4], Br: [20, 11, -3], I: [-6, 11, -1], OH: [49, 10, -6], OR: [49, 7, -5], OPh: [54, 6, -5], OAc: [51, 6, -3], NH2: [29, 11, -5], NHR: [37, 8, -4], NR2: [42, 6, -3],
  NHAc: [28, 7, -4], NR3: [26, 8, -5], NArom: [30, 6, -3], NO2: [63, 4, -4], SH: [11, 12, -4], SR: [20, 7, -3], SO: [31, 5, -3], SO2: [30, 5, -3], CN: [4, 3, -3], vinyl: [20, 6, -0.5], alkynyl: [4.5, 5.5, -3.5],
  Ph: [23, 9, -2], CHO: [30, -0.5, -2.5], COR: [30, 1, -2], COOH: [21, 3, -2], COOR: [20, 3, -2], CON: [22, 0, -3], COCl: [33, 3, -2], Si: [-3, 2, -1], P: [10, 5, -2], B: [0, 0, 0], X: [0, 0, 0] };
/* alkene carbon increments [on the carbon the group is attached to, on the other alkene carbon] */
const C_ALKENE = { Ph: [12.5, -11.0], vinyl: [13.6, -7.0], alkynyl: [-6.0, 6.0], CHO: [13.1, 12.7], COR: [4.5, 6.0], COOH: [5.0, 9.8], COOR: [6.3, 7.0], CON: [5.0, 7.0], COCl: [8.1, 14.0], CN: [-15.1, 14.2],
  OH: [28.8, -39.5], OR: [28.8, -39.5], OPh: [28.8, -39.5], OAc: [18.4, -26.7], NH2: [28.0, -32.0], NHR: [28.0, -32.0], NR2: [28.0, -32.0], NHAc: [23.0, -20.0], NArom: [23.0, -20.0], NO2: [22.3, -0.9],
  F: [24.9, -34.3], Cl: [2.0, -6.1], Br: [-8.0, -0.9], I: [-38.1, 7.0], SH: [18.0, -16.0], SR: [18.0, -16.0], SO2: [14.5, 7.1], SO: [14.5, 7.1], Si: [16.0, -7.0] };

const DAMP_H = [1, 0.8, 0.7, 0.6], DAMP_C = [1, 0.85, 0.7, 0.6];
function damped(vals, damp) { const v = vals.slice().sort((a, b) => Math.abs(b) - Math.abs(a)); return v.reduce((s, x, i) => s + x * damp[Math.min(i, damp.length - 1)], 0); }

/* ------------------------------------------------------------------ ring geometry helpers for aromatic positions */
function ringOrder(r) { return r.atoms; }                     /* findRings returns atoms in cycle order */
function ringDist(r, a, b) { const n = r.size, i = r.atoms.indexOf(a), j = r.atoms.indexOf(b), d = Math.abs(i - j); return Math.min(d, n - d); }
function homeRing(G, i) { let best = -1; G.atoms[i].rings.forEach(k => { if (best < 0 || G.rings[k].size < G.rings[best].size) best = k; }); return best; }
const isPyN = (G, i) => { const a = G.atoms[i]; return a.el === 'N' && a.arom && a.nH === 0 && a.chg <= 0 && a.nb.some(x => x.bd.o === 2 || x.bd.ar) && a.nb.length === 2; };
const isPyrrN = (G, i) => { const a = G.atoms[i]; return a.arom && ((a.el === 'N' && !isPyN(G, i)) || a.el === 'O' || a.el === 'S' || a.el === 'Se'); };

/* Substituent increments for an aromatic ring position (works for 1H and 13C). */
function arSCS(G, ring, c, tab, nIdx) {
  let s = 0; const parts = [];
  for (const r of ring.atoms) {
    const d = ringDist(ring, c, r);
    for (const { j } of G.atoms[r].nb) {
      if (ring.atoms.includes(j)) continue;
      const key = subKey(G, r, j), row = tab[key] || (G.atoms[j].el === 'C' ? tab.Et : null); if (!row) continue;
      let pos = d, v;
      if (nIdx === 4) { v = row[Math.min(d, 3)]; }                 /* 13C: index 0 ipso, 1 ortho, 2 meta, 3 para */
      else { if (d === 0) continue; v = row[Math.min(d, 3) - 1]; }
      if (ring.size === 5 && d === 2 && nIdx === 3) v = row[1];
      /* in a fused ring system the ring-fusion atom's partner ring acts as an aryl substituent with a smaller ipso effect */
      if (nIdx === 4 && d === 0 && key === 'Ph' && G.atoms[r].rings.length > 1) v = 5.2;
      s += v; parts.push([key, d, v]);
    }
  }
  return s;
}

function nEffect(G, ring, c, which) {                        /* pyridine-type N in a six-membered ring: alpha, beta, gamma */
  const A = which === 'H' ? [1.33, -0.02, 0.37] : [21.4, -4.7, 7.4], damp2 = which === 'H' ? 0.7 : 0.35;
  let total = 0, alphaSeen = 0;
  for (const r of ring.atoms) if (r !== c && isPyN(G, r)) {
    const d = ringDist(ring, c, r) - 1;
    if (d === 0) { total += A[0] * (alphaSeen ? damp2 : 1); alphaSeen++; } else total += A[d];
  }
  if (isPyN(G, c)) return 0;
  return total;
}

/* Base value for a carbon or hydrogen on a five-membered heteroaromatic ring. */
function fiveRingBase(G, ring, c, which) {
  const ringHet = ring.atoms.filter(i => G.atoms[i].el !== 'C');
  const pos = ring.atoms.indexOf(c), n = ring.size;
  const nb = [ring.atoms[(pos + 1) % n], ring.atoms[(pos + n - 1) % n]], far = [ring.atoms[(pos + 2) % n], ring.atoms[(pos + n - 2) % n]];
  const adjN = nb.filter(i => isPyN(G, i)).length, adjX = nb.filter(i => isPyrrN(G, i)).length;
  const X = ringHet.map(i => G.atoms[i]).find(a => isPyrrN(G, G.atoms.indexOf(a))); const xe = X ? (X.el === 'N' ? 'N' : X.el) : 'N';
  const H = which === 'H';
  if (ringHet.length === 1) {
    const e = G.atoms[ringHet[0]].el, alpha = nb.some(i => G.atoms[i].el !== 'C');
    const T = H ? { O: [7.40, 6.35], S: [7.30, 7.08], N: [6.68, 6.22] } : { O: [142.8, 109.8], S: [125.6, 127.3], N: [118.0, 108.0] };
    return T[e === 'O' || e === 'S' ? e : 'N'][alpha ? 0 : 1];
  }
  if (adjN && adjX) return (H ? { S: 8.75, O: 7.9, N: 7.7 } : { S: 153, O: 151, N: 136 })[xe];
  if (adjN === 2) return H ? 8.0 : 140;
  if (adjN) { const x2 = far.some(i => isPyrrN(G, i)); return H ? (x2 ? 7.55 : { S: 7.85, O: 7.1, N: 7.1 }[xe]) : (x2 ? 139 : { S: 143, O: 125, N: 122 }[xe]); }
  if (adjX) { const n2 = far.some(i => isPyN(G, i)); return H ? (n2 ? 7.5 : { S: 7.3, O: 7.7, N: 7.1 }[xe]) : (n2 ? 129.5 : { S: 120, O: 138, N: 122 }[xe]); }
  return H ? 6.3 : 105;
}

/* ------------------------------------------------------------------ 1H prediction */
function alphaIncH(G, c, s) {
  const k = subKey(G, c, s);
  if ((k === 'NH2' || k === 'NHR' || k === 'NR2') && hasArNbr(G, s, c)) return 1.9;
  return H_ALPHA[k] || 0;
}

function predictH(G, opt) {
  const out = [], solvent = (opt && opt.solvent) || 'CDCl3';
  for (let i = 0; i < G.n; i++) {
    const a = G.atoms[i]; if (a.nH <= 0) continue;
    const base = { atom: i, cls: a.cls, el: a.el };
    if (a.el !== 'C') { for (let k = 0; k < a.nH; k++) out.push(Object.assign({}, base, exchangeH(G, i, solvent))); continue; }
    if (a.arom) { const r = G.rings[homeRing(G, i)]; const v = aromaticH(G, i, r); out.push(Object.assign({}, base, { kind: 'arom', shift: v, label: 'ArH', env: 'ring' })); continue; }
    const o = dblO(G, i);
    if (o >= 0) { /* aldehyde */
      const oth = a.nb.filter(x => x.j !== o)[0];
      let v = 9.75; if (oth) { const e = G.atoms[oth.j]; v = e.arom ? 10.0 : (e.el === 'N' || e.el === 'O') ? 8.05 : isAlkeneC(G, oth.j) ? 9.6 : 9.75; }
      out.push(Object.assign({}, base, { kind: 'ald', shift: v, label: 'CHO' })); continue;
    }
    const t = tripleNbr(G, i);
    if (t >= 0) { out.push(Object.assign({}, base, { kind: 'alkyne', shift: alkyneH(G, i, t), label: '≡CH' })); continue; }
    if (a.nb.some(x => x.bd.o === 2)) {
      const p = a.nb.find(x => x.bd.o === 2).j, pe = G.atoms[p].el;
      if (pe !== 'C') { out.push(Object.assign({}, base, { kind: 'imine', shift: 8.0, label: 'HC=N' })); continue; }
      for (const h of vinylH(G, i, p)) out.push(Object.assign({}, base, h)); continue;
    }
    if (a.chg !== 0) continue;
    const v = sp3H(G, i); for (let k = 0; k < a.nH; k++) out.push(Object.assign({}, base, { kind: 'sp3', shift: v, label: 'CH' + (a.nH > 1 ? a.nH : '') }));
  }
  /* merge equivalent protons into signals */
  const sig = new Map();
  for (const h of out) {
    const key = h.cls + '|' + h.kind + '|' + (h.side === undefined ? '' : h.shiftKey || '') + '|' + h.shift.toFixed(2);
    let s = sig.get(key); if (!s) { s = { shift: h.shift, kind: h.kind, label: h.label, atoms: [], hs: [], cls: h.cls, exch: h.kind === 'exch' }; sig.set(key, s); }
    s.hs.push(h); if (!s.atoms.includes(h.atom)) s.atoms.push(h.atom);
  }
  const list = [...sig.values()];
  list.forEach((s, k) => { s.id = k; s.n = s.hs.length; s.hs.forEach(h => h.sig = k); });
  couple(G, out, list);
  list.sort((x, y) => y.shift - x.shift);
  return { protons: out, signals: list };
}

function exchangeH(G, i, solvent) {
  const a = G.atoms[i], nb = a.nb.map(x => x.j), dm = solvent === 'DMSO-d6';
  let v = 2.5, env = '';
  if (a.el === 'O') {
    const c = nb[0];
    if (c === undefined) { v = 4.8; env = 'water-like'; }
    else if (isCarbonyl(G, c) || (G.atoms[c].el === 'S') || (G.atoms[c].el === 'P')) { v = dm ? 12.2 : 11.0; env = 'acid O–H'; }
    else if (G.atoms[c].arom) { v = dm ? 9.3 : 5.5; env = 'phenol O–H'; }
    else if (G.atoms[c].el === 'C' && hasPi(G, c) && isAlkeneC(G, c)) { v = dm ? 8.0 : 5.5; env = 'enol O–H'; }
    else { v = dm ? 4.5 : 2.5; env = 'alcohol O–H'; }
  } else if (a.el === 'N') {
    const c = nb.find(j => true);
    if (a.arom) { v = dm ? 11.0 : 8.0; env = 'ring N–H'; }
    else if (a.chg > 0) { v = dm ? 8.5 : 8.0; env = 'N⁺–H'; }
    else if (c !== undefined && nb.some(j => isCarbonyl(G, j))) { const prim = a.nH === 2; v = (prim ? 6.5 : 7.8) + (dm ? 0.8 : 0); env = 'amide N–H'; if (nb.some(j => isCarbonyl(G, j) && G.atoms[j].nb.some(y => G.atoms[y.j].el === 'O' && y.bd.o === 1) || G.atoms[j].nb.filter(y => G.atoms[y.j].el === 'N').length > 1)) { v = 5.5 + (dm ? 0.8 : 0); env = 'carbamate / urea N–H'; } }
    else if (nb.some(j => G.atoms[j].el === 'S' && G.atoms[j].nb.some(y => y.bd.o === 2))) { v = 7.8; env = 'sulfonamide N–H'; }
    else if (nb.some(j => G.atoms[j].arom)) { v = (a.nH === 2 ? 3.6 : 3.8) + (dm ? 1.2 : 0); env = 'aryl amine N–H'; }
    else { v = (a.nH === 2 ? 1.3 : 1.5) + (dm ? 0.5 : 0); env = 'amine N–H'; }
  } else if (a.el === 'S') {
    v = nb.some(j => G.atoms[j].arom) ? 3.4 : 1.4; env = 'thiol S–H';
  } else { v = 2.0; env = a.el + '–H'; }
  return { kind: 'exch', shift: v, label: a.el + 'H', env, side: undefined };
}

function sp3H(G, c) {
  const a = G.atoms[c], n = a.nH, ring = smallestRing(G, c);
  let base = n === 3 ? 0.90 : n === 2 ? 1.27 : 1.55;
  if (ring && n === 2) base = RING_BASE_CH2[ring] || 1.5; else if (ring && n === 1) base = RING_BASE_CH[ring] || 1.6;
  const al = []; let beta = 0;
  for (const { j } of a.nb) {
    if (isSp3C(G, j)) {
      for (const { j: e } of G.atoms[j].nb) if (e !== c) { const k = subKey(G, j, e); beta += H_BETA[k] || 0; }
    } else al.push(alphaIncH(G, c, j));
  }
  return base + damped(al, DAMP_H) + beta;
}

function alkyneH(G, i, t) {
  const o = G.atoms[t].nb.filter(x => x.j !== i)[0]; if (!o) return 1.8;
  const k = subKey(G, t, o.j);
  return k === 'Ph' ? 3.05 : /^CO|^CHO|vinyl/.test(k) ? 2.9 : /^CH2[ONX]/.test(k) ? 2.4 : 1.9;
}

function sideOf(G, i, p, s) { const A = G.atoms[i], P = G.atoms[p], S = G.atoms[s]; const z = (P.x - A.x) * (S.y - A.y) - (P.y - A.y) * (S.x - A.x); return z >= 0 ? 1 : -1; }
function vinylH(G, i, p) {
  const a = G.atoms[i], heavy = a.nb.filter(x => x.j !== p), pSubs = G.atoms[p].nb.filter(x => x.j !== i).map(x => x.j);
  const Hs = [];
  if (heavy.length === 0) Hs.push({ side: 1 }, { side: -1 });
  else if (heavy.length === 1) Hs.push({ side: -sideOf(G, i, p, heavy[0].j) });
  const gem = heavy.length ? subKey(G, i, heavy[0].j) : null;
  return Hs.map(h => {
    let v = 5.25; if (gem) v += (H_VINYL[gem] || ALK_ALKYL)[0];
    for (const u of pSubs) { const k = subKey(G, p, u), row = H_VINYL[k] || (G.atoms[u].el === 'C' ? ALK_ALKYL : [0, 0, 0]); v += sideOf(G, i, p, u) === h.side ? row[1] : row[2]; }
    return { kind: 'vinyl', shift: v, label: '=CH', side: h.side, shiftKey: String(h.side), partner: p };
  });
}

function aromaticH(G, c, ring) {
  let v;
  if (ring.size === 6) v = 7.27 + nEffect(G, ring, c, 'H'); else v = fiveRingBase(G, ring, c, 'H');
  return v + arSCS(G, ring, c, AR_H, 3);
}

/* ------------------------------------------------------------------ couplings and multiplets */
const J_VIC = 7.0, J_VIC_ALLYL = 6.5, J_ALD = 2.5, J_CIS = 10.5, J_TRANS = 16.8, J_GEM_VINYL = 1.5;
function ringJ(G, i, j) {                                   /* aromatic H on i and j, which are ring neighbours */
  const k = homeRing(G, i), r = G.rings[k], six = r.size === 6;
  if (six) { const ni = r.atoms.some(x => isPyN(G, x)), al = x => G.atoms[x].nb.some(y => isPyN(G, y.j)); return ni && (al(i) !== al(j)) ? 4.8 : 7.8; }
  const het = r.atoms.map(x => G.atoms[x]).find(a => a.el !== 'C'); const e = het ? het.el : 'N';
  const alphaPair = [i, j].some(x => G.atoms[x].nb.some(y => G.atoms[y.j].el !== 'C' && r.atoms.includes(y.j)));
  return alphaPair ? (e === 'O' ? 1.8 : e === 'S' ? 5.0 : 2.7) : 3.4;
}
function couple(G, protons, signals) {
  const byAtom = new Map(); protons.forEach(h => { if (!byAtom.has(h.atom)) byAtom.set(h.atom, []); byAtom.get(h.atom).push(h); });
  for (const sig of signals) {
    sig.J = []; sig.coupled = [];
    if (sig.exch) continue;
    const rep = sig.hs[0], i = rep.atom, a = G.atoms[i], partners = [];
    const add = (h, J) => { if (h.sig !== sig.id && !h.exch && h.kind !== 'exch') partners.push({ sig: h.sig, J }); };
    /* geminal vinyl */
    if (rep.kind === 'vinyl') for (const h of byAtom.get(i)) if (h !== rep) add(h, J_GEM_VINYL);
    for (const { j } of a.nb) {
      const hs = byAtom.get(j); if (!hs) continue;
      for (const h of hs) {
        if (h.kind === 'exch') continue;
        let J = 0;
        if (rep.kind === 'arom' && h.kind === 'arom') J = ringJ(G, i, j);
        else if (rep.kind === 'vinyl' && h.kind === 'vinyl' && G.atoms[i].nb.some(x => x.j === j && x.bd.o === 2)) J = rep.side === h.side ? J_CIS : J_TRANS;
        else if (rep.kind === 'ald' || h.kind === 'ald') J = (rep.kind === 'sp3' || h.kind === 'sp3') ? J_ALD : 0;
        else if (rep.kind === 'sp3' && h.kind === 'sp3') J = J_VIC;
        else if ((rep.kind === 'sp3' && h.kind === 'vinyl') || (rep.kind === 'vinyl' && h.kind === 'sp3')) J = J_VIC_ALLYL;
        if (J > 0) add(h, J);
      }
    }
    /* meta coupling on six-membered aromatic rings */
    if (rep.kind === 'arom') {
      const r = G.rings[homeRing(G, i)];
      if (r.size === 6) for (const x of r.atoms) if (ringDist(r, i, x) === 2) { const hs = byAtom.get(x); if (hs && hs[0].kind === 'arom') for (const h of hs) add(h, 1.8); }
    }
    /* count partner protons per (signal, J) */
    const m = new Map(); partners.forEach(p => { const k = p.sig + '|' + p.J; const e = m.get(k) || { sig: p.sig, J: p.J, n: 0 }; e.n++; m.set(k, e); });
    sig.coupled = [...m.values()];
  }
  for (const sig of signals) {
    if (sig.exch) { sig.mult = 'br s'; sig.lines = [{ d: 0, w: 1 }]; sig.Jlist = []; continue; }
    /* merge couplings with nearly equal J */
    const groups = []; sig.coupled.slice().sort((a, b) => b.J - a.J).forEach(c => { const g = groups.find(x => Math.abs(x.J - c.J) < 0.8); if (g) { g.J = (g.J * g.n + c.J * c.n) / (g.n + c.n); g.n += c.n; } else groups.push({ J: c.J, n: c.n }); });
    sig.Jlist = groups; sig.lines = multiplet(groups); sig.mult = label(groups, sig.lines.length);
  }
}
function binom(n, k) { let r = 1; for (let i = 1; i <= k; i++) r = r * (n - k + i) / i; return r; }
function multiplet(groups) {
  let lines = [{ d: 0, w: 1 }];
  for (const g of groups) {
    const nl = []; for (const l of lines) for (let k = 0; k <= g.n; k++) nl.push({ d: l.d + (k - g.n / 2) * g.J, w: l.w * binom(g.n, k) });
    lines = nl;
  }
  lines.sort((a, b) => a.d - b.d); const out = [];
  for (const l of lines) { const p = out[out.length - 1]; if (p && Math.abs(p.d - l.d) < 0.05) p.w += l.w; else out.push({ d: l.d, w: l.w }); }
  const tot = out.reduce((s, l) => s + l.w, 0); out.forEach(l => l.w /= tot);
  return out;
}
const NAME1 = ['s', 'd', 't', 'q', 'quint', 'sext', 'sept'], LET = ['', 'd', 't', 'q'];
function label(groups, nlines) {
  if (!groups.length) return 's';
  if (groups.length === 1) return groups[0].n <= 6 ? NAME1[groups[0].n] : 'm';
  if (groups.length === 2 && groups.every(g => g.n <= 3)) return LET[groups[0].n] + LET[groups[1].n];
  if (groups.length === 3 && groups.every(g => g.n === 1)) return 'ddd';
  return 'm';
}

/* ------------------------------------------------------------------ 13C prediction */
const CORR = { 1: { 3: -1.1, 4: -3.4 }, 2: { 3: -2.5, 4: -7.5 }, 3: { 2: -3.7, 3: -9.5 }, 4: { 1: -1.5, 2: -8.4 } };
const RING_CORR = { 3: -18.7, 4: -2.9, 5: -8.8, 6: -5.2, 7: -1.0, 8: -3.0 };
function cSp3(G, c) {
  const dist = new Map([[c, 0]]), q = [c]; let nA = 0, nB = 0, nG = 0, nD = 0; const al = [], be = [], ga = [];
  while (q.length) {
    const u = q.shift(), du = dist.get(u); if (du >= 4) continue;
    for (const { j: v } of G.atoms[u].nb) {
      if (dist.has(v)) continue; dist.set(v, du + 1);
      if (isSp3C(G, v)) { q.push(v); const d = du + 1; if (d === 1) nA++; else if (d === 2) nB++; else if (d === 3) nG++; else nD++; }
      else if (du + 1 <= 3) { const k = subKey(G, u, v), row = C_SUB[k]; if (row) [al, be, ga][du].push(row[du]); }
    }
  }
  let v = -2.3 + 9.1 * nA + 9.4 * nB - 2.5 * nG + 0.3 * nD + damped(al, DAMP_C) + be.reduce((s, x) => s + x, 0) + ga.reduce((s, x) => s + x, 0);
  const nC = G.atoms[c].nb.filter(x => isSp3C(G, x.j)).length;
  if (CORR[nC]) for (const { j } of G.atoms[c].nb) if (isSp3C(G, j)) { const dj = G.atoms[j].nb.filter(x => isSp3C(G, x.j)).length; v += CORR[nC][dj] || 0; }
  const k = smallestRing(G, c);
  if (k) { const r = G.rings.filter(r => !r.arom && r.atoms.includes(c) && r.size === k)[0]; if (r && r.atoms.every(x => isSp3C(G, x))) v += RING_CORR[k] || 0; }
  return v;
}

function cAlkene(G, c) {
  const a = G.atoms[c], pn = a.nb.find(x => x.bd.o === 2), p = pn.j; let v = 123.3;
  const side = (i, other, tabIdx) => {
    for (const { j: s } of G.atoms[i].nb) {
      if (s === other) continue;
      if (isSp3C(G, s)) {
        const A = tabIdx === 0 ? [10.6, 7.2, -1.5] : [-7.9, -1.8, 1.5]; v += A[0];
        for (const { j: b } of G.atoms[s].nb) if (b !== i && isSp3C(G, b)) { v += A[1]; for (const { j: g } of G.atoms[b].nb) if (g !== s && isSp3C(G, g)) v += A[2]; }
      } else { const k = subKey(G, i, s), row = C_ALKENE[k]; if (row) v += row[tabIdx]; }
    }
  };
  side(c, p, 0); side(p, c, 1);
  return v;
}

function cArom(G, c) {
  const r = G.rings[homeRing(G, c)];
  let base; if (r.size === 6) base = 128.5 + nEffect(G, r, c, 'C'); else base = fiveRingBase(G, r, c, 'C');
  return base + arSCS(G, r, c, AR_C, 4);
}

function cCarbonyl(G, c) {
  const a = G.atoms[c], o = dblO(G, c), rest = a.nb.filter(x => x.j !== o);
  const cl = x => { const e = G.atoms[x.j]; if (e.el === 'C') return e.arom ? 'r' : (isAlkeneC(G, x.j) || tripleNbr(G, x.j) >= 0 ? 'v' : 'C'); if (e.el === 'O') return e.chg < 0 || e.nH > 0 ? 'OH' : e.nb.some(y => y.j !== c && isCarbonyl(G, y.j)) ? 'OAcyl' : 'O'; return e.el; };
  const k = rest.map(cl).sort(), ring = G.atoms[c].nb.length && smallestRing(G, c);
  const pair = { 'C,C': 208, 'C,r': 198, 'C,v': 199, 'r,r': 196, 'r,v': 191, 'v,v': 189 };
  const C3 = { C: 0, r: 1, v: 1 };
  const nH = a.nH;
  if (nH > 0) { if (!k.length) return 200; const t = k[0]; return t === 'r' ? 192 : t === 'v' ? 193 : t === 'N' ? 162 : t === 'O' ? 161 : 201; }
  let v = null;
  if (k.length === 2 && k.every(t => /^[Crv]$/.test(t))) { v = pair[k.join(',')]; if (ring) v += { 4: 0, 5: 10, 6: 3, 7: 7 }[ring] || 0; return v; }
  if (k.length === 1) {
    const t = k[0]; if (t === 'Cl') return 170; if (t === 'S') return 198;
  }
  const hasC = k.find(t => /^[Crv]$/.test(t)), other = k.find(t => !/^[Crv]$/.test(t));
  if (hasC && other) {
    const sh = hasC === 'C' ? 0 : hasC === 'r' ? -1 : -1;
    if (other === 'OH') return (hasC === 'C' ? 178 : 172);
    if (other === 'O') { v = hasC === 'C' ? 172 : 166; if (ring === 5) v += 5; return v; }
    if (other === 'OAcyl') return 168;
    if (other === 'N') { v = hasC === 'C' ? 173 : hasC === 'r' ? 168 : 166; if (ring === 5) v += 2; return v; }
    if (other === 'Cl') return 170;
    if (other === 'S') return 198;
  }
  if (k.length === 2) { const s = k.join(','); if (/O/.test(s) && /N/.test(s)) return 156; if (s === 'N,N') return 160; if (/O.*O|OH.*O/.test(s)) return 155; if (/Cl/.test(s)) return 150; }
  return 175;
}

function cSp(G, c) {
  const a = G.atoms[c], t = tripleNbr(G, c), te = G.atoms[t].el;
  if (te === 'N') return 118.5;
  const own = a.nb.filter(x => x.j !== t)[0], oth = G.atoms[t].nb.filter(x => x.j !== c)[0];
  const kind = (n, from) => { if (!n) return 'H'; const k = subKey(G, from, n.j); return k === 'Ph' ? 'aryl' : /^(CO|CHO|CON)/.test(k) ? 'CO' : 'alkyl'; };
  const A = { H: -12, alkyl: 0, aryl: 1, CO: -3 }, B = { H: 4, alkyl: 0, aryl: 8.7, CO: -4 };
  const k1 = G.atoms[c].nH > 0 ? 'H' : kind(own, c), k2 = G.atoms[t].nH > 0 ? 'H' : kind(oth, t);
  return 80.5 + A[k1] + B[k2];
}

function predictC(G, opt) {
  const out = [];
  for (let i = 0; i < G.n; i++) {
    const a = G.atoms[i]; if (a.el !== 'C') continue;
    let v, type;
    if (a.arom) { v = cArom(G, i); }
    else if (isCarbonyl(G, i)) v = cCarbonyl(G, i);
    else if (tripleNbr(G, i) >= 0) v = cSp(G, i);
    else if (a.nb.some(x => x.bd.o === 2)) { const p = a.nb.find(x => x.bd.o === 2).j; v = G.atoms[p].el === 'C' ? (a.nb.filter(x => x.bd.o === 2).length > 1 ? 100 : cAlkene(G, i)) : 158; }
    else if (a.chg === 0) v = cSp3(G, i);
    else v = 100;
    type = a.nH === 3 ? 'CH3' : a.nH === 2 ? 'CH2' : a.nH === 1 ? 'CH' : 'C';
    out.push({ atom: i, cls: a.cls, shift: v, type, nH: a.nH });
  }
  const sig = new Map();
  for (const c of out) {
    const key = c.cls + '|' + c.shift.toFixed(1); let s = sig.get(key);
    if (!s) { s = { shift: c.shift, type: c.type, atoms: [], n: 0, nH: c.nH }; sig.set(key, s); }
    s.atoms.push(c.atom); s.n++;
  }
  const list = [...sig.values()].sort((x, y) => y.shift - x.shift); list.forEach((s, k) => s.id = k);
  return { carbons: out, signals: list };
}

/* ------------------------------------------------------------------ top level */
function formula(G) {
  const cnt = {}; G.atoms.forEach(a => { cnt[a.el] = (cnt[a.el] || 0) + 1; cnt.H = (cnt.H || 0) + a.nH; });
  const order = ['C', 'H'].concat(Object.keys(cnt).filter(k => k !== 'C' && k !== 'H').sort());
  return order.filter(k => cnt[k]).map(k => k + (cnt[k] > 1 ? cnt[k] : '')).join('');
}
function predict(molfileText, opt) {
  const mol = parseMolfile(molfileText), G = build(mol);
  if (G.unknown.length) throw new Error('Not supported yet: ' + [...new Set(G.unknown)].join(', '));
  const H = predictH(G, opt), C = predictC(G, opt);
  const nHtot = G.atoms.reduce((s, a) => s + a.nH, 0), nHobs = H.signals.filter(s => !s.exch).reduce((s, x) => s + x.n, 0);
  return { G, formula: formula(G), nHtotal: nHtot, nHobserved: nHobs, h: H, c: C };
}

/* ------------------------------------------------------------------ spectrum synthesis */
/* lines: [{ppm, w}] where w is the area of the line. Lorentzian shape, width at half height lw (Hz) at spectrometer frequency mhz. */
function synth(lines, o) {
  const { xmin, xmax, N, lw, mhz } = o, g = (lw / 2) / mhz, y = new Float64Array(N), dx = (xmax - xmin) / (N - 1), reach = Math.max(40 * g, 0.2);
  for (const l of lines) {
    const i0 = Math.max(0, Math.floor((l.ppm - reach - xmin) / dx)), i1 = Math.min(N - 1, Math.ceil((l.ppm + reach - xmin) / dx));
    for (let i = i0; i <= i1; i++) { const x = xmin + i * dx - l.ppm; y[i] += l.w * (g / Math.PI) / (x * x + g * g); }
  }
  return y;
}
/* 1H line list: every multiplet line, area proportional to the number of protons. */
function hLines(sigs, mhz, opt) {
  const out = [], showEx = !opt || opt.exchangeable !== false;
  for (const s of sigs) {
    if (s.exch && !showEx) continue;
    const bw = s.exch ? 4 : 0;
    for (const l of s.lines) out.push({ ppm: s.shift + l.d / mhz, w: s.n * l.w, sig: s.id, broad: bw });
  }
  return out;
}
function cLines(sigs, mode) {
  const out = [];
  for (const s of sigs) {
    if (mode === 'dept') { if (s.type === 'C') continue; out.push({ ppm: s.shift, w: s.n * (s.type === 'CH2' ? -1 : 1), sig: s.id }); }
    else out.push({ ppm: s.shift, w: s.n * (s.type === 'C' ? 0.45 : s.type === 'CH' ? 0.9 : 1), sig: s.id });
  }
  return out;
}

const API = { parseMolfile, build, predict, synth, hLines, cLines, formula, version: '1' };
if (typeof module !== 'undefined' && module.exports) module.exports = API; else root.NMR = API;
})(typeof self !== 'undefined' ? self : this);
