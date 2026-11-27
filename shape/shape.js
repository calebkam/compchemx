/* ============================ shape.js ============================
   CCX Shape: Lewis structures, VSEPR geometry and molecular shape. Pure functions, no DOM, no network.

   Input : a V2000 molfile (from CCX Draw).
   Output: Lewis structure (valence electrons, lone pairs, formal charges, resonance forms, octet check),
           VSEPR analysis per atom (AXnEm, electron geometry, molecular shape, hybridisation, ideal and measured angles),
           a 3D model, electron-domain directions for drawing clouds, and a bond-dipole polarity estimate.

   3D geometry
     * Acyclic molecules, and anything with a hypervalent or lone-pair centre, are built here by VSEPR: each atom's
       electron domains (bonds, lone pairs, unpaired electrons) sit on the ideal polyhedron and are relaxed under
       weighted repulsion, so lone pairs squeeze bond angles (H2O about 104.5, NH3 about 107) and axial / equatorial
       positions come out right (SF4 seesaw, ClF3 T-shaped, XeF4 square planar).
     * Rings are handed to the force-field builder (Chem.build3D from dock/chem.js) and the electron domains are then
       derived from that geometry.
   Bond lengths come from covalent radii and bond order, so they are good to about 0.05 A. VSEPR is a model for main-group
   molecules. It does not describe transition-metal complexes.                                                        */
(function (root) {
'use strict';
const PI = Math.PI, R2D = 180 / PI, D2R = PI / 180;

/* ve valence electrons, en Pauling electronegativity, r single-bond covalent radius (A), per period, val allowed valences */
const EL = {
  H: { ve: 1, en: 2.20, r: .32, per: 1 }, He: { ve: 2, en: 0, r: .46, per: 1 },
  Be: { ve: 2, en: 1.57, r: .96, per: 2 }, B: { ve: 3, en: 2.04, r: .85, per: 2 }, C: { ve: 4, en: 2.55, r: .75, per: 2 },
  N: { ve: 5, en: 3.04, r: .71, per: 2 }, O: { ve: 6, en: 3.44, r: .63, per: 2 }, F: { ve: 7, en: 3.98, r: .64, per: 2 }, Ne: { ve: 8, en: 0, r: .67, per: 2 },
  Al: { ve: 3, en: 1.61, r: 1.21, per: 3 }, Si: { ve: 4, en: 1.90, r: 1.16, per: 3 }, P: { ve: 5, en: 2.19, r: 1.11, per: 3 }, S: { ve: 6, en: 2.58, r: 1.03, per: 3 },
  Cl: { ve: 7, en: 3.16, r: .99, per: 3 }, Ar: { ve: 8, en: 0, r: .96, per: 3 },
  Ge: { ve: 4, en: 2.01, r: 1.20, per: 4 }, As: { ve: 5, en: 2.18, r: 1.19, per: 4 }, Se: { ve: 6, en: 2.55, r: 1.16, per: 4 }, Br: { ve: 7, en: 2.96, r: 1.14, per: 4 }, Kr: { ve: 8, en: 3.00, r: 1.17, per: 4 },
  Sn: { ve: 4, en: 1.96, r: 1.39, per: 5 }, Sb: { ve: 5, en: 2.05, r: 1.38, per: 5 }, Te: { ve: 6, en: 2.10, r: 1.35, per: 5 }, I: { ve: 7, en: 2.66, r: 1.33, per: 5 }, Xe: { ve: 8, en: 2.60, r: 1.40, per: 5 },
};
const VAL_BY_GROUP = { 0: [0], 1: [1], 2: [2], 3: [3], 4: [4], 5: [3, 5], 6: [2, 4, 6], 7: [1, 3, 5, 7], 8: [0, 2, 4, 6] };
const known = el => !!EL[el];

/* ------------------------------------------------------------------ molfile */
function parseMolfile(text) {
  const ln = String(text).replace(/\r/g, '').split('\n');
  if (!/V2000/.test(ln[3] || '')) throw new Error('Expected a V2000 molfile.');
  const na = parseInt(ln[3].substr(0, 3)), nb = parseInt(ln[3].substr(3, 3));
  if (!(na > 0)) throw new Error('The drawing has no atoms.');
  const CH = { 1: 3, 2: 2, 3: 1, 5: -1, 6: -2, 7: -3 }, atoms = [], bonds = [];
  for (let i = 0; i < na; i++) { const l = ln[4 + i]; atoms.push({ el: l.substr(31, 3).trim(), x: parseFloat(l.substr(0, 10)) || 0, y: parseFloat(l.substr(10, 10)) || 0, chg: CH[parseInt(l.substr(36, 3))] || 0 }); }
  for (let i = 0; i < nb; i++) { const l = ln[4 + na + i]; bonds.push({ a: parseInt(l.substr(0, 3)) - 1, b: parseInt(l.substr(3, 3)) - 1, o: Math.min(3, Math.max(1, parseInt(l.substr(6, 3)) || 1)) }); }
  let saw = false;
  for (const l of ln) if (/^M  CHG/.test(l)) { if (!saw) { atoms.forEach(a => a.chg = 0); saw = true; } const n = parseInt(l.substr(6, 3)); for (let i = 0; i < n; i++) { const k = parseInt(l.substr(10 + i * 8, 4)) - 1, v = parseInt(l.substr(14 + i * 8, 4)); if (atoms[k]) atoms[k].chg = v; } }
  return { atoms, bonds };
}

/* ------------------------------------------------------------------ graph, hydrogens, lone pairs */
function validValences(el, q) {
  const e = EL[el], g = e.ve - q; let list = VAL_BY_GROUP[Math.max(0, Math.min(8, g))] || [0];
  if (el === 'H') return q === 0 ? [1] : [0];
  if (e.per === 2) list = list.filter(v => v <= 4); if (!list.length) list = [0];
  return list;
}
function build(mol) {
  const atoms = mol.atoms.map((a, i) => ({ el: a.el === 'D' || a.el === 'T' ? 'H' : a.el, x: a.x, y: a.y, chg: a.chg || 0, src: i, imp: false, nb: [] }));
  const bonds = [];
  for (const b of mol.bonds) { if (b.a === b.b) continue; const bd = { a: b.a, b: b.b, o: b.o, i: bonds.length }; bonds.push(bd); atoms[b.a].nb.push({ j: b.b, bd }); atoms[b.b].nb.push({ j: b.a, bd }); }
  const warnings = [], n0 = atoms.length;
  for (const a of atoms) if (!known(a.el)) throw new Error('Not supported yet: ' + a.el + '. This tool covers main-group elements (H, B to F, Al to Cl, Ge to Br, Sn to Xe).');
  for (let i = 0; i < n0; i++) {
    const a = atoms[i]; let bo = 0; for (const { bd } of a.nb) bo += bd.o;
    const v = validValences(a.el, a.chg).find(x => x >= bo);
    let nH = 0; if (v === undefined) warnings.push(a.el + (i + 1) + ' has more bonds than a normal valence allows.'); else nH = v - bo;
    for (let k = 0; k < nH; k++) {
      const h = { el: 'H', x: a.x, y: a.y, chg: 0, src: -1, imp: true, nb: [] }, bd = { a: i, b: atoms.length, o: 1, i: bonds.length };
      atoms.push(h); bonds.push(bd); a.nb.push({ j: atoms.length - 1, bd }); h.nb.push({ j: i, bd });
    }
  }
  const G = { atoms, bonds, n: atoms.length, nHeavy: n0, warnings };
  for (const a of atoms) {
    const e = EL[a.el]; let bo = 0; for (const { bd } of a.nb) bo += bd.o;
    const free = e.ve - a.chg - bo;           /* electrons left over after bonding */
    a.bo = bo; a.lp = Math.max(0, Math.floor(free / 2)); a.rad = free % 2 !== 0 && free > 0 ? 1 : 0;
    if (free < 0) warnings.push(a.el + ' has more bonding electrons than it has valence electrons.');
  }
  G.cycles = cycles(G); G.hasRing = G.bonds.length - G.n + components(G) > 0;
  return G;
}
function components(G) { const seen = new Array(G.n).fill(false); let c = 0; for (let i = 0; i < G.n; i++) if (!seen[i]) { c++; const q = [i]; seen[i] = true; while (q.length) { const u = q.pop(); for (const { j } of G.atoms[u].nb) if (!seen[j]) { seen[j] = true; q.push(j); } } } return c; }
function cycles(G, maxSize = 8) {
  const out = [], seen = new Set(), n = G.n;
  for (let s = 0; s < n && out.length < 300; s++) {
    if (G.atoms[s].el === 'H') continue; const path = [s], on = new Array(n).fill(false); on[s] = true;
    (function dfs(u) {
      for (const { j } of G.atoms[u].nb) {
        if (G.atoms[j].el === 'H') continue;
        if (j === s && path.length >= 3) { const k = path.slice().sort((a, b) => a - b).join(','); if (!seen.has(k)) { seen.add(k); out.push(path.slice()); } continue; }
        if (j <= s || on[j] || path.length >= maxSize) continue; path.push(j); on[j] = true; dfs(j); path.pop(); on[j] = false;
      }
    })(s);
  }
  return out;
}

/* ------------------------------------------------------------------ resonance */
function stateOf(G) { return { bo: G.bonds.map(b => b.o), lp: G.atoms.map(a => a.lp), q: G.atoms.map(a => a.chg) }; }
const keyOf = s => s.bo.join('') + '|' + s.lp.join('') + '|' + s.q.join('');
function bondOf(G, a, b) { for (const { j, bd } of G.atoms[a].nb) if (j === b) return bd.i; return -1; }
function moves(G, s) {
  const out = [], at = G.atoms;
  const mk = ch => { const t = { bo: s.bo.slice(), lp: s.lp.slice(), q: s.q.slice() }; ch(t); return t; };
  const bsum = (t, i) => at[i].nb.reduce((x, { bd }) => x + t.bo[bd.i], 0);
  for (let X = 0; X < G.n; X++) {
    if (at[X].el === 'H' || s.lp[X] < 1) continue;                      /* M1: lone pair on X pushes into X-Y, Y=Z becomes Y-Z with a lone pair on Z */
    for (const { j: Y, bd: b1 } of at[X].nb) { if (s.bo[b1.i] >= 3 || at[Y].el === 'H') continue;
      for (const { j: Z, bd: b2 } of at[Y].nb) { if (Z === X || s.bo[b2.i] < 2 || at[Z].el === 'H') continue;
        out.push(mk(t => { t.bo[b1.i]++; t.bo[b2.i]--; t.lp[X]--; t.lp[Z]++; t.q[X]++; t.q[Z]--; })); } }
  }
  for (let Z = 0; Z < G.n; Z++) {                                          /* M2: electron-poor Z (cation, B) takes a pi bond from Y=X */
    const e = EL[at[Z].el]; if (e.per !== 2 || at[Z].el === 'H') continue;
    if (2 * bsum(s, Z) + 2 * s.lp[Z] >= 8) continue; if (!(s.q[Z] > 0 || at[Z].el === 'B')) continue;
    for (const { j: Y, bd: b2 } of at[Z].nb) { if (s.bo[b2.i] >= 3) continue;
      for (const { j: X, bd: b1 } of at[Y].nb) { if (X === Z || s.bo[b1.i] < 2 || at[X].el === 'H') continue;
        out.push(mk(t => { t.bo[b2.i]++; t.bo[b1.i]--; t.q[Z]--; t.q[X]++; })); } }
  }
  for (const cyc of G.cycles) {                                             /* M3: flip an alternating ring (Kekule forms) */
    if (cyc.length % 2) continue; const bi = cyc.map((x, k) => bondOf(G, x, cyc[(k + 1) % cyc.length])); if (bi.some(b => b < 0)) continue;
    const o = bi.map(b => s.bo[b]); if (!(o.every((v, k) => v === (k % 2 ? 1 : 2)) || o.every((v, k) => v === (k % 2 ? 2 : 1)))) continue;
    out.push(mk(t => { bi.forEach(b => t.bo[b] = 3 - t.bo[b]); }));
  }
  return out;
}
function penalty(G, s) {
  let p = 0;
  G.atoms.forEach((a, i) => {
    const e = EL[a.el];
    if (e.per === 2 && a.el !== 'H' && a.el !== 'B' && a.el !== 'Be') { const el = 2 * a.nb.reduce((x, { bd }) => x + s.bo[bd.i], 0) + 2 * s.lp[i]; if (el < 8) p += 10; }
    p += 2 * Math.abs(s.q[i]);
    if (s.q[i] < 0) p += (3.98 - e.en) * 0.4 * -s.q[i]; if (s.q[i] > 0) p += e.en * 0.2 * s.q[i];
  });
  return p;
}
function resonance(G, limit = 60) {
  const s0 = stateOf(G), seen = new Map([[keyOf(s0), s0]]), queue = [s0];
  while (queue.length && seen.size < limit) { const s = queue.shift(); for (const t of moves(G, s)) { const k = keyOf(t); if (!seen.has(k)) { seen.set(k, t); queue.push(t); if (seen.size >= limit) break; } } }
  const forms = [...seen.values()].map((s, i) => ({ state: s, p: penalty(G, s), drawn: i === 0 }));
  const pmin = Math.min(...forms.map(f => f.p)); forms.forEach(f => f.major = f.p <= pmin + 0.01);
  forms.sort((a, b) => (b.drawn - a.drawn) || (a.p - b.p));
  const major = forms.filter(f => f.major), hyb = { bo: G.bonds.map((b, k) => major.reduce((x, f) => x + f.state.bo[k], 0) / major.length), q: G.atoms.map((a, k) => major.reduce((x, f) => x + f.state.q[k], 0) / major.length) };
  const minLP = G.atoms.map((a, k) => Math.min(...forms.map(f => f.state.lp[k])));
  return { forms, major, hybrid: hyb, minLP, count: forms.length, majorCount: major.length };
}

/* ------------------------------------------------------------------ Lewis summary */
function lewis(G, R) {
  const total = G.atoms.reduce((s, a) => s + EL[a.el].ve, 0) - G.atoms.reduce((s, a) => s + a.chg, 0);
  const net = G.atoms.reduce((s, a) => s + a.chg, 0), notes = [];
  const bondE = G.bonds.reduce((s, b) => s + 2 * b.o, 0), lpE = G.atoms.reduce((s, a) => s + 2 * a.lp, 0), unp = G.atoms.reduce((s, a) => s + a.rad, 0);
  const atoms = G.atoms.map((a, i) => {
    const e = EL[a.el], bo = a.bo, el = 2 * bo + 2 * a.lp + a.rad;
    let octet = 'ok';
    if (a.el === 'H') octet = el === 2 ? 'ok' : el < 2 ? 'incomplete' : 'over'; else if (a.el === 'He') octet = 'ok';
    else if (el === 8) octet = 'ok'; else if (el > 8) octet = e.per >= 3 ? 'expanded' : 'over'; else octet = (a.el === 'B' || a.el === 'Be' || a.el === 'Al') && a.chg === 0 ? 'deficient (normal for ' + a.el + ')' : 'incomplete';
    return { i, el: a.el, q: a.chg, lp: a.lp, rad: a.rad, bondingE: 2 * bo, electrons: el, octet, fc: a.chg };
  });
  if (unp) notes.push('Odd number of electrons: ' + unp + ' unpaired electron' + (unp > 1 ? 's' : '') + ' (radical).');
  atoms.filter(a => a.octet === 'expanded').forEach(a => notes.push(a.el + (a.i + 1) + ' has an expanded octet (' + a.electrons + ' electrons), which period 3 and later atoms allow.'));
  atoms.filter(a => a.octet === 'over').forEach(a => notes.push(a.el + (a.i + 1) + ' exceeds an octet, which a period 2 atom cannot do.'));
  atoms.filter(a => a.octet === 'incomplete' && a.el !== 'H').forEach(a => notes.push(a.el + (a.i + 1) + ' has an incomplete octet (' + a.electrons + ' electrons).'));
  if (R && R.count > 1) notes.push(R.count + ' resonance form' + (R.count > 1 ? 's' : '') + ' found (' + R.majorCount + ' major).');
  return { totalValence: total, pairs: total / 2, bondingElectrons: bondE, lonePairElectrons: lpE, unpaired: unp, netCharge: net, atoms, notes, formula: formula(G) };
}
function formula(G) {
  const c = {}; G.atoms.forEach(a => c[a.el] = (c[a.el] || 0) + 1);
  const keys = Object.keys(c), order = c.C ? ['C', 'H'].concat(keys.filter(k => k !== 'C' && k !== 'H').sort()) : keys.sort((a, b) => a === 'H' ? -1 : b === 'H' ? 1 : a < b ? -1 : 1);
  const net = G.atoms.reduce((s, a) => s + a.chg, 0);
  return order.filter(k => c[k]).map(k => k + (c[k] > 1 ? c[k] : '')).join('') + (net ? (Math.abs(net) > 1 ? Math.abs(net) : '') + (net > 0 ? '+' : '−') : '');
}

/* ------------------------------------------------------------------ VSEPR tables */
const SHAPES = {
  '10': ['terminal atom', 'terminal atom', ''], '20': ['linear', 'linear', '180°'], '21': ['trigonal planar', 'bent', 'a little under 120°'], '22': ['tetrahedral', 'bent', 'a little under 109.5° (water 104.5°)'],
  '23': ['trigonal bipyramidal', 'linear', '180°'], '30': ['trigonal planar', 'trigonal planar', '120°'], '31': ['tetrahedral', 'trigonal pyramidal', 'a little under 109.5° (ammonia 107°)'],
  '32': ['trigonal bipyramidal', 'T-shaped', 'about 87.5° and 175°'], '33': ['octahedral', 'T-shaped', '90°, 180°'], '40': ['tetrahedral', 'tetrahedral', '109.5°'], '41': ['trigonal bipyramidal', 'seesaw', 'about 102° (equatorial), 173° (axial)'],
  '42': ['octahedral', 'square planar', '90°, 180°'], '50': ['trigonal bipyramidal', 'trigonal bipyramidal', '90°, 120°, 180°'], '51': ['octahedral', 'square pyramidal', 'a little under 90°'],
  '52': ['pentagonal bipyramidal', 'pentagonal planar', '72°, 144°'], '60': ['octahedral', 'octahedral', '90°, 180°'], '61': ['pentagonal bipyramidal', 'pentagonal pyramidal', '72°, 90°'],
  '70': ['pentagonal bipyramidal', 'pentagonal bipyramidal', '72°, 90°, 180°'],
};
const HYB = { 1: '', 2: 'sp', 3: 'sp²', 4: 'sp³', 5: 'sp³d', 6: 'sp³d²', 7: 'sp³d³' };
const s3 = Math.sqrt(3) / 2, t4 = 1 / Math.sqrt(3);
const TEMPLATE = {
  1: [[0, 0, 1]], 2: [[0, 0, 1], [0, 0, -1]], 3: [[1, 0, 0], [-.5, s3, 0], [-.5, -s3, 0]],
  4: [[t4, t4, t4], [t4, -t4, -t4], [-t4, t4, -t4], [-t4, -t4, t4]],
  5: [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-.5, s3, 0], [-.5, -s3, 0]],
  6: [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]],
  7: [[0, 0, 1], [0, 0, -1]].concat([0, 1, 2, 3, 4].map(k => [Math.cos(k * 2 * PI / 5), Math.sin(k * 2 * PI / 5), 0])),
};
const LP_SLOTS = { 2: { 1: [1] }, 3: { 1: [2], 2: [1, 2] }, 4: { 1: [3], 2: [2, 3], 3: [1, 2, 3] }, 5: { 1: [2], 2: [2, 3], 3: [2, 3, 4] }, 6: { 1: [1], 2: [0, 1], 3: [0, 1, 5], 4: [0, 1, 4, 5] }, 7: { 1: [0], 2: [0, 1] } };
const EQ_FIRST = { 5: [2, 3, 4, 0, 1] };                 /* trigonal bipyramid: put the bulky ligands equatorial */
const W_LP = 1.22, W_RAD = 0.85, W_DBL = 1.16, W_TRI = 1.3;

const v3 = { add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]], sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], mul: (a, k) => [a[0] * k, a[1] * k, a[2] * k],
  dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2], cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
  len: a => Math.hypot(a[0], a[1], a[2]), nrm: a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; } };
function rotAxis(v, axis, ang) { const k = v3.nrm(axis), c = Math.cos(ang), s = Math.sin(ang), d = v3.dot(k, v), x = v3.cross(k, v); return [v[0] * c + x[0] * s + k[0] * d * (1 - c), v[1] * c + x[1] * s + k[1] * d * (1 - c), v[2] * c + x[2] * s + k[2] * d * (1 - c)]; }
function alignRot(a, b) {                                  /* rotation taking unit a to unit b, as a function */
  const c = v3.dot(a, b), ax = v3.cross(a, b), l = v3.len(ax);
  if (l < 1e-9) { if (c > 0) return v => v; let p = Math.abs(a[0]) < .9 ? [1, 0, 0] : [0, 1, 0]; p = v3.nrm(v3.cross(a, p)); return v => rotAxis(v, p, PI); }
  const ang = Math.atan2(l, c); return v => rotAxis(v, ax, ang);
}
function relax(P, W, iters = 400) {
  P = P.map(p => v3.nrm(p)); const n = P.length;
  for (let it = 0; it < iters; it++) {
    const g = P.map(() => [0, 0, 0]);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const d = v3.sub(P[i], P[j]), r2 = v3.dot(d, d) || 1e-9, f = -2 * W[i] * W[j] / (r2 * r2);
      for (let k = 0; k < 3; k++) { g[i][k] += f * d[k]; g[j][k] -= f * d[k]; }
    }
    for (let i = 0; i < n; i++) { const rad = v3.dot(g[i], P[i]); const t = [g[i][0] - rad * P[i][0], g[i][1] - rad * P[i][1], g[i][2] - rad * P[i][2]]; P[i] = v3.nrm([P[i][0] - 0.03 * t[0], P[i][1] - 0.03 * t[1], P[i][2] - 0.03 * t[2]]); }
  }
  return P;
}

/* Domains of one atom: neighbours (sigma domains weighted by bond order), effective lone pairs, unpaired electrons. */
function domainsOf(G, i, R, hybOrder) {
  const a = G.atoms[i], lpEff = R ? Math.max(0, Math.min(a.lp, R.minLP[i])) : a.lp;
  const nb = a.nb.map(({ j, bd }) => ({ j, o: Math.round(hybOrder ? hybOrder[bd.i] : bd.o) || 1, bi: bd.i }));
  return { nb, X: nb.length, E: lpEff, lp: a.lp, rad: a.rad, SN: nb.length + lpEff + a.rad, deloc: a.lp - lpEff };
}
function template(dom, nbInfo) {
  const SN = dom.SN, E = dom.E + dom.rad; if (SN < 1 || SN > 7) return null;
  const P0 = TEMPLATE[SN], lpSlots = (LP_SLOTS[SN] && LP_SLOTS[SN][E]) || (E ? null : []);
  if (lpSlots === null) return null;
  const ligSlots = P0.map((_, k) => k).filter(k => !lpSlots.includes(k));
  const order = (EQ_FIRST[SN] || ligSlots).filter(k => ligSlots.includes(k));
  const lig = nbInfo.map((n, idx) => ({ idx, w: n.o === 3 ? W_TRI : n.o === 2 ? W_DBL : 1, en: EL[n.el].en })).sort((p, q) => (q.w - p.w) || (p.en - q.en));
  const slotOf = {}; lig.forEach((l, k) => slotOf[l.idx] = order[k]);
  const W = P0.map(() => W_LP); lig.forEach(l => W[slotOf[l.idx]] = l.w); lpSlots.forEach((s, k) => W[s] = k < dom.E ? W_LP : W_RAD);
  const P = relax(P0, W);
  return { P, slotOf, lpSlots: lpSlots.slice(), SN };
}
function bondLength(G, a, b, order, hypervalent) {
  const A = EL[G.atoms[a].el], B = EL[G.atoms[b].el]; let f = 1;
  if (order > 1) f = order <= 2 ? 1 - 0.11 * (order - 1) : 0.89 - 0.09 * (order - 2);
  if (hypervalent) f *= 0.95;
  return (A.r + B.r) * f;
}

/* ------------------------------------------------------------------ 3D builder (VSEPR, acyclic) */
function dihedral(p0, p1, p2, p3) {
  const b0 = v3.sub(p0, p1), b1 = v3.nrm(v3.sub(p2, p1)), b2 = v3.sub(p3, p2);
  const v = v3.sub(b0, v3.mul(b1, v3.dot(b0, b1))), w = v3.sub(b2, v3.mul(b1, v3.dot(b2, b1)));
  return Math.atan2(v3.dot(v3.cross(b1, v), w), v3.dot(v, w));
}
function buildVSEPR(G, R) {
  const n = G.n, X = new Array(n).fill(null), parent = new Array(n).fill(-1), done = new Array(n).fill(false);
  const hybOrder = R ? R.hybrid.bo : null, doms = [];
  for (let i = 0; i < n; i++) doms.push(domainsOf(G, i, R, hybOrder));
  let root = 0; for (let i = 1; i < n; i++) if (doms[i].SN > doms[root].SN || (doms[i].SN === doms[root].SN && G.atoms[i].el !== 'H' && G.atoms[root].el === 'H')) root = i;
  const lpWorld = Array.from({ length: n }, () => []), radWorld = Array.from({ length: n }, () => []);
  const sp2like = i => doms[i].SN <= 3 && (G.atoms[i].nb.some(({ bd }) => (hybOrder ? hybOrder[bd.i] : bd.o) > 1.2) || doms[i].deloc > 0);
  const hyper = i => doms[i].SN >= 5;
  X[root] = [0, 0, 0]; const queue = [root]; done[root] = true;
  while (queue.length) {
    const a = queue.shift(), d = doms[a], nbInfo = d.nb.map(nb => ({ ...nb, el: G.atoms[nb.j].el }));
    const T = template(d, nbInfo);
    if (!T) { if (d.SN > 7) throw new Error('This atom has more than seven electron domains.'); if (d.SN > 0) throw new Error('No VSEPR template for ' + G.atoms[a].el + (a + 1) + ' (' + d.X + ' bonds, ' + (d.E + d.rad) + ' non-bonding).'); continue; }
    /* orientation: align the parent's slot onto the bond back to the parent */
    let rot = v => v, par = parent[a];
    if (par >= 0) { const pi = d.nb.findIndex(x => x.j === par), tgt = v3.nrm(v3.sub(X[par], X[a])); rot = alignRot(T.P[T.slotOf[pi]], tgt); }
    const children = d.nb.map((nb, idx) => ({ ...nb, idx })).filter(c => c.j !== par && !done[c.j]);
    const placeAt = (phi) => {
      const ax = par >= 0 ? v3.nrm(v3.sub(X[par], X[a])) : [0, 0, 1], r = v => par >= 0 ? rotAxis(rot(v), ax, phi) : v;
      return { dirs: T.P.map(r), r };
    };
    let best = null, cands = [];
    if (par < 0) cands = [0]; else {
      /* planar candidates when both ends are flat pi systems, else staggered by clash search */
      const base = placeAt(0), g = G.atoms[par].nb.map(x => x.j).find(j => j !== a && X[j]);
      const first = children[0];
      const planarWanted = first && g !== undefined && sp2like(a) && sp2like(par) && (G.atoms[a].nb.find(x => x.j === par).bd.o > 1 || (hybOrder && hybOrder[G.atoms[a].nb.find(x => x.j === par).bd.i] > 1.2) || doms[a].deloc > 0 || doms[par].deloc > 0);
      if (planarWanted) {
        const p1 = X[par], p2 = X[a], p3 = v3.add(X[a], v3.mul(base.dirs[T.slotOf[first.idx]], 1)), t0 = dihedral(X[g], p1, p2, p3);
        cands = [-t0, PI - t0];
      } else for (let k = 0; k < 12; k++) cands.push(k * PI / 6);
    }
    let bestScore = -1;
    for (const phi of cands) {
      const { dirs } = placeAt(phi); let min = 1e9;
      const pts = [];
      children.forEach(c => { const L = bondLength(G, a, c.j, hybOrder ? hybOrder[c.bi] : c.o, hyper(a) || hyper(c.j)); pts.push(v3.add(X[a], v3.mul(dirs[T.slotOf[c.idx]], L))); });
      T.lpSlots.forEach(s => pts.push(v3.add(X[a], dirs[s])));
      for (let k = 0; k < n; k++) { if (!X[k] || k === a || k === par) continue; for (const p of pts) min = Math.min(min, v3.len(v3.sub(p, X[k]))); }
      for (let u = 0; u < pts.length; u++) for (let w = u + 1; w < pts.length; w++) { /* keep own virtual points apart from nothing; ignored */ }
      if (par >= 0) for (const lw of lpWorld[par]) for (const p of pts) min = Math.min(min, v3.len(v3.sub(p, v3.add(X[par], lw))) * 0.8);
      const sc = min + (planarWanted0(cands, phi) ? 0 : 0);
      if (sc > bestScore + 1e-9) { bestScore = sc; best = { phi, dirs }; }
    }
    function planarWanted0() { return 0; }
    const { dirs } = best;
    children.forEach(c => { const L = bondLength(G, a, c.j, hybOrder ? hybOrder[c.bi] : c.o, hyper(a) || hyper(c.j)); X[c.j] = v3.add(X[a], v3.mul(dirs[T.slotOf[c.idx]], L)); parent[c.j] = a; done[c.j] = true; queue.push(c.j); });
    T.lpSlots.forEach((s, k) => { (k < d.E ? lpWorld[a] : radWorld[a]).push(dirs[s]); });
    /* delocalised lone pairs are shown as pi: remember the true pair count separately */
  }
  /* lone pairs that were delocalised: re-add as perpendicular pi-type pair (drawn as a p lobe) */
  for (let i = 0; i < n; i++) if (!X[i]) throw new Error('The molecule has separate pieces. Draw one molecule at a time.');
  return { X, lpWorld, radWorld, doms };
}

/* Lone-pair directions around an atom whose neighbours are already placed (used for ring systems). */
function placeLonePairs(center, nbDirs, k, weightsNb) {
  if (k <= 0) return [];
  const nb = nbDirs.map(d => v3.nrm(d)); let bestE = 1e18, best = null;
  for (let seed = 0; seed < 24; seed++) {
    let lp = []; for (let m = 0; m < k; m++) { const a = (seed * 2.399963 + m * 2.094395) % (2 * PI), b = Math.acos(1 - 2 * (((seed * 0.618 + m * 0.382 + 0.07) % 1))); lp.push([Math.sin(b) * Math.cos(a), Math.sin(b) * Math.sin(a), Math.cos(b)]); }
    const W = nb.map((_, i) => (weightsNb && weightsNb[i]) || 1).concat(lp.map(() => W_LP)), all = () => nb.concat(lp);
    for (let it = 0; it < 250; it++) {
      const P = all(), g = lp.map(() => [0, 0, 0]);
      for (let i = 0; i < k; i++) for (let j = 0; j < P.length; j++) { if (j === nb.length + i) continue; const d = v3.sub(P[nb.length + i], P[j]), r2 = v3.dot(d, d) || 1e-9, f = -2 * W[nb.length + i] * W[j] / (r2 * r2); for (let c = 0; c < 3; c++) g[i][c] += f * d[c]; }
      lp = lp.map((p, i) => { const rad = v3.dot(g[i], p), t = [g[i][0] - rad * p[0], g[i][1] - rad * p[1], g[i][2] - rad * p[2]]; return v3.nrm([p[0] - 0.03 * t[0], p[1] - 0.03 * t[1], p[2] - 0.03 * t[2]]); });
    }
    const P = all(); let E = 0; for (let i = 0; i < P.length; i++) for (let j = i + 1; j < P.length; j++) { const d = v3.sub(P[i], P[j]); E += W[i] * W[j] / (v3.dot(d, d) || 1e-9); }
    if (E < bestE - 1e-9) { bestE = E; best = lp; }
  }
  return best;
}

/* ------------------------------------------------------------------ analysis from coordinates */
function angleDeg(P, i, a, b) { const u = v3.nrm(v3.sub(P[a], P[i])), w = v3.nrm(v3.sub(P[b], P[i])); return Math.acos(Math.max(-1, Math.min(1, v3.dot(u, w)))) * R2D; }
function analyse(G, R, P, lpWorld) {
  return G.atoms.map((a, i) => {
    const d = domainsOf(G, i, R, R ? R.hybrid.bo : null), key = '' + d.X + (d.E + d.rad), tab = SHAPES[key] || (d.SN > 7 ? ['—', '—', ''] : ['—', 'unusual arrangement', '']);
    const nbs = a.nb.map(x => x.j), angs = []; for (let p = 0; p < nbs.length; p++) for (let q = p + 1; q < nbs.length; q++) angs.push({ i: nbs[p], k: nbs[q], deg: angleDeg(P, i, nbs[p], nbs[q]) });
    const lengths = a.nb.map(({ j, bd }) => ({ j, d: v3.len(v3.sub(P[i], P[j])), o: R ? R.hybrid.bo[bd.i] : bd.o }));
    const terminal = d.X <= 1;
    return { i, el: a.el, X: d.X, E: d.E, rad: d.rad, SN: d.SN, lonePairsDrawn: a.lp, delocalised: d.deloc,
      axe: terminal ? '' : 'AX' + (d.X > 1 ? d.X : '') + (d.E + d.rad ? 'E' + (d.E + d.rad > 1 ? d.E + d.rad : '') : ''), electronGeometry: terminal ? '' : tab[0], shape: terminal ? 'terminal atom' : tab[1], ideal: tab[2],
      hybridisation: terminal ? '' : HYB[d.SN] || '', angles: angs, lengths, center: !terminal };
  });
}
function polarity(G, R, P) {
  const net = [0, 0, 0], bonds = [];
  G.bonds.forEach(b => {
    const A = EL[G.atoms[b.a].el], B = EL[G.atoms[b.b].el], d = (B.en - A.en); if (Math.abs(d) < 1e-6) { bonds.push({ a: b.a, b: b.b, d: 0 }); return; }
    const toward = d > 0 ? b.b : b.a, from = d > 0 ? b.a : b.b, u = v3.nrm(v3.sub(P[toward], P[from])), m = Math.abs(d);
    for (let k = 0; k < 3; k++) net[k] += u[k] * m; bonds.push({ a: from, b: toward, d: m });
  });
  const tot = G.atoms.reduce((x, a) => x + a.chg, 0);
  if (R && tot === 0) { const c = P.reduce((x, p) => v3.add(x, p), [0, 0, 0]).map(v => v / P.length); G.atoms.forEach((a, i) => { const q = R.hybrid.q[i]; if (Math.abs(q) > 1e-6) for (let k = 0; k < 3; k++) net[k] += 0.5 * q * (P[i][k] - c[k]); }); }
  const mag = v3.len(net), nonzero = bonds.filter(b => b.d > 0.05).length;
  return { vector: net, magnitude: mag, polar: mag > 0.2, bonds, hasPolarBonds: nonzero > 0, estimateDebye: mag };
}

/* ------------------------------------------------------------------ 2D Lewis layout (hydrogens and lone-pair dots placed in free angles) */
function layout2D(G, state) {
  const A = G.atoms, st = state || stateOf(G), pos = A.map(a => ({ x: a.x, y: -a.y }));
  const lens = []; G.bonds.forEach(b => { if (!A[b.a].imp && !A[b.b].imp) lens.push(Math.hypot(pos[b.a].x - pos[b.b].x, pos[b.a].y - pos[b.b].y)); });
  lens.sort((x, y) => x - y); const bl = lens.length ? lens[Math.floor(lens.length / 2)] : 1.5 > 0 ? (lens[0] || 1) : 1;
  const L = bl || 1;
  const ang = (i, j) => Math.atan2(pos[j].y - pos[i].y, pos[j].x - pos[i].x);
  const free = i => {
    const a = A[i], placed = a.nb.filter(x => !A[x.j].imp || pos[x.j].placed).map(x => ang(i, x.j)).sort((p, q) => p - q);
    return placed;
  };
  const gapsFill = (i, count) => {
    const out = []; let angles = A[i].nb.filter(x => !A[x.j].imp || A[x.j].placed).map(x => ang(i, x.j));
    for (let k = 0; k < count; k++) {
      const s = angles.slice().sort((p, q) => p - q);
      if (!s.length) { out.push(-PI / 2 + k * 2 * PI / count); angles.push(out[out.length - 1]); continue; }
      let bestGap = -1, mid = 0;
      for (let m = 0; m < s.length; m++) { const a0 = s[m], a1 = m + 1 < s.length ? s[m + 1] : s[0] + 2 * PI, gap = a1 - a0; if (gap > bestGap + 1e-9) { bestGap = gap; mid = a0 + gap / 2; } }
      out.push(mid); angles.push(mid);
    }
    return out;
  };
  for (let i = 0; i < G.nHeavy; i++) {
    const hs = A[i].nb.filter(x => A[x.j].imp).map(x => x.j);
    if (!hs.length) continue; const dirs = gapsFill(i, hs.length);
    hs.forEach((h, k) => { pos[h] = { x: pos[i].x + Math.cos(dirs[k]) * L * 0.85, y: pos[i].y + Math.sin(dirs[k]) * L * 0.85 }; A[h].placed = true; });
  }
  const dots = [];
  for (let i = 0; i < G.n; i++) {
    const lp = st.lp[i], rad = A[i].rad; if (!lp && !rad) continue;
    const dirs = gapsFill(i, lp + rad);
    dirs.forEach((d, k) => {
      const c = [pos[i].x + Math.cos(d) * L * 0.36, pos[i].y + Math.sin(d) * L * 0.36];
      if (k < lp) { const px = -Math.sin(d) * L * 0.07, py = Math.cos(d) * L * 0.07; dots.push({ x: c[0] + px, y: c[1] + py, atom: i }, { x: c[0] - px, y: c[1] - py, atom: i }); } else dots.push({ x: c[0], y: c[1], atom: i });
    });
  }
  return { pos, dots, bl: L, bonds: G.bonds.map(b => ({ a: b.a, b: b.b, o: st.bo[b.i], i: b.i })), q: st.q };
}

/* ------------------------------------------------------------------ top level */
function analyze(molfileText, opt) {
  opt = opt || {};
  const mol = parseMolfile(molfileText), G = build(mol), R = resonance(G), L = lewis(G, R);
  let X = null, lpWorld = null, radWorld = null, method = 'vsepr', warn = G.warnings.slice();
  let needVSEPR = true;
  if (G.hasRing) {
    if (!opt.chem) throw new Error('Ring systems need the force-field builder (dock/chem.js), which is not loaded.');
    const coords = opt.chem.build3D(molfileText, { nconf: 6, seed: 11 });
    /* the force-field model has its own hydrogens: match them to ours by heavy-atom parent order */
    const ff = coords.atoms, heavyFF = []; ff.forEach((a, i) => { if (a.el !== 'H') heavyFF.push(i); });
    if (heavyFF.length !== G.nHeavy) throw new Error('The force-field model has a different atom count than the Lewis structure.');
    X = new Array(G.n); const map = new Array(G.n).fill(-1); for (let i = 0; i < G.nHeavy; i++) { map[i] = heavyFF[i]; X[i] = [ff[heavyFF[i]].x, ff[heavyFF[i]].y, ff[heavyFF[i]].z]; }
    const used = new Set(); for (let i = 0; i < G.nHeavy; i++) {
      const hs = G.atoms[i].nb.filter(x => G.atoms[x.j].imp).map(x => x.j);
      const cand = coords.bonds.map(b => b.a === map[i] ? b.b : b.b === map[i] ? b.a : -1).filter(j => j >= 0 && ff[j].el === 'H' && !used.has(j));
      hs.forEach((h, k) => { const j = cand[k]; used.add(j); X[h] = [ff[j].x, ff[j].y, ff[j].z]; });
    }
    for (let i = 0; i < G.n; i++) if (!X[i]) X[i] = [0, 0, 0];
    method = 'force field (ring system)'; lpWorld = G.atoms.map(() => []); radWorld = G.atoms.map(() => []);
    for (let i = 0; i < G.n; i++) { const d = domainsOf(G, i, R, R.hybrid.bo); if (d.E + d.rad > 0 && d.X > 0) { const nbd = G.atoms[i].nb.map(x => v3.nrm(v3.sub(X[x.j], X[i]))); const dirs = placeLonePairs(X[i], nbd, d.E + d.rad, G.atoms[i].nb.map(x => R.hybrid.bo[x.bd.i] > 1.5 ? W_DBL : 1)); dirs.forEach((u, k) => (k < d.E ? lpWorld[i] : radWorld[i]).push(u)); } }
    needVSEPR = false;
  }
  let doms;
  if (needVSEPR) { const b = buildVSEPR(G, R); X = b.X; lpWorld = b.lpWorld; radWorld = b.radWorld; doms = b.doms; }
  /* terminal atoms (one neighbour) with lone pairs: place them around the bond, staggered */
  for (let i = 0; i < G.n; i++) if (needVSEPR && G.atoms[i].nb.length === 1 && !lpWorld[i].length && G.atoms[i].lp > 0) {
    const d = domainsOf(G, i, R, R.hybrid.bo), nbd = [v3.nrm(v3.sub(X[G.atoms[i].nb[0].j], X[i]))];
    if (d.E + d.rad > 0) placeLonePairs(X[i], nbd, d.E + d.rad, [R.hybrid.bo[G.atoms[i].nb[0].bd.i] > 1.5 ? W_DBL : 1]).forEach((u, k) => (k < d.E ? lpWorld[i] : radWorld[i]).push(u));
  }
  const cen = X.reduce((s, p) => v3.add(s, p), [0, 0, 0]).map(v => v / G.n); X = X.map(p => v3.sub(p, cen));
  const atoms = analyse(G, R, X, lpWorld), pol = polarity(G, R, X);
  const model = { atoms: G.atoms.map((a, i) => ({ el: a.el, x: X[i][0], y: X[i][1], z: X[i][2], q: a.chg, imp: a.imp, lp: lpWorld[i], rad: radWorld[i], lpCount: a.lp, en: EL[a.el].en })), bonds: G.bonds.map(b => ({ a: b.a, b: b.b, o: b.o, order: R.hybrid.bo[b.i] })) };
  return { G, resonance: R, lewis: L, atoms, polarity: pol, model, method, warnings: warn, layout: layout2D(G, R.forms[0].state) };
}

const API = { parseMolfile, build, resonance, lewis, analyze, layout2D, placeLonePairs, relax, SHAPES, EL, formula, stateOf, version: '1' };
if (typeof module !== 'undefined' && module.exports) module.exports = API; else root.Shape = API;
})(typeof self !== 'undefined' ? self : this);
