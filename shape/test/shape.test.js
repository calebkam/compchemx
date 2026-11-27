/* Shape tests: node shape/test/shape.test.js */
const S = require('../shape.js'), C = require('../../dock/chem.js'), F = require('./molfiles.json');
let bad = 0, n = 0;
const ok = (c, m) => { n++; if (!c) { bad++; console.log('FAIL', m); } };
const near = (a, b, t, m) => ok(Math.abs(a - b) <= t, m + ' got ' + (+a).toFixed(1) + ' want ' + b + '±' + t);
const run = k => S.analyze(F[k], { chem: C });
const ctr = (r, el) => r.atoms.find(a => a.center && a.el === el);
const ang = (a, deg, t = 3) => a.angles.some(g => Math.abs(g.deg - deg) <= t);
const SHAPE = { H2O: ['O', 'AX2E2', 'bent'], NH3: ['N', 'AX3E', 'trigonal pyramidal'], CH4: ['C', 'AX4', 'tetrahedral'], BF3: ['B', 'AX3', 'trigonal planar'], SF4: ['S', 'AX4E', 'seesaw'], ClF3: ['Cl', 'AX3E2', 'T-shaped'], XeF2: ['Xe', 'AX2E3', 'linear'], XeF4: ['Xe', 'AX4E2', 'square planar'], SF6: ['S', 'AX6', 'octahedral'], PCl5: ['P', 'AX5', 'trigonal bipyramidal'], 'ICl4-': ['I', 'AX4E2', 'square planar'], CO2: ['C', 'AX2', 'linear'], SO2: ['S', 'AX2E', 'bent'], 'NH4+': ['N', 'AX4', 'tetrahedral'] };
for (const k in SHAPE) { const a = ctr(run(k), SHAPE[k][0]); ok(a && a.axe === SHAPE[k][1] && a.shape === SHAPE[k][2], k + ' VSEPR ' + (a && a.axe + ' ' + a.shape)); }
near(ctr(run('H2O'), 'O').angles[0].deg, 104.5, 2, 'H2O angle'); near(ctr(run('NH3'), 'N').angles[0].deg, 107, 2, 'NH3 angle');
ok(ang(ctr(run('SF6'), 'S'), 90, 1) && ang(ctr(run('SF6'), 'S'), 180, 1), 'SF6 90/180');
ok(ang(ctr(run('XeF4'), 'Xe'), 90, 1) && ang(ctr(run('XeF4'), 'Xe'), 180, 1), 'XeF4 90/180');
ok(ang(ctr(run('SF4'), 'S'), 173, 8) && ang(ctr(run('SF4'), 'S'), 87, 5), 'SF4 axial/eq angles');
ok(ang(ctr(run('BF3'), 'B'), 120, 1), 'BF3 120');
ok(ang(ctr(run('ClF3'), 'Cl'), 87, 5), 'ClF3 ~87');
/* Lewis */
const L = k => run(k).lewis;
ok(L('H2O').totalValence === 8 && L('CO3').totalValence === 24 && L('SF6').totalValence === 48 && L('XeF4').totalValence === 36 && L('NH4+').totalValence === 8, 'valence electron totals');
ok(run('XeF4').G.atoms[0].lp + run('XeF4').G.atoms.find(a => a.el === 'Xe').lp >= 2 && run('XeF4').G.atoms.find(a => a.el === 'Xe').lp === 2, 'XeF4 Xe has 2 lone pairs');
ok(run('CO3').lewis.netCharge === -2 && run('CO3').resonance.majorCount === 3, 'carbonate: 3 major forms');
ok(run('benzene').resonance.majorCount === 2, 'benzene 2 Kekule forms');
ok(run('nitromethane').resonance.majorCount === 2, 'nitro 2 forms');
ok(run('SF6').lewis.atoms.find(a => a.el === 'S').octet === 'expanded', 'SF6 expanded octet');
ok(run('BF3').lewis.atoms.find(a => a.el === 'B').octet.startsWith('deficient'), 'BF3 deficient');
/* polarity */
for (const k of ['CCl4', 'CO2', 'BF3', 'SF6', 'XeF4', 'CH4', 'XeF2', 'benzene', 'ethene', 'PCl5']) ok(!run(k).polarity.polar, k + ' nonpolar');
for (const k of ['H2O', 'NH3', 'CHCl3', 'SO2', 'SF4', 'ClF3', 'ozone', 'ethanol']) ok(run(k).polarity.polar, k + ' polar');
/* geometry sanity: no two atoms closer than 0.6 A */
for (const k in F) { const m = run(k).model; let min = 9; for (let i = 0; i < m.atoms.length; i++) for (let j = i + 1; j < m.atoms.length; j++) min = Math.min(min, Math.hypot(m.atoms[i].x - m.atoms[j].x, m.atoms[i].y - m.atoms[j].y, m.atoms[i].z - m.atoms[j].z)); ok(min > 0.6, k + ' atoms overlap ' + min.toFixed(2)); }
console.log(n - bad + '/' + n + ' passed'); process.exit(bad ? 1 : 0);
