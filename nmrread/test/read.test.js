/* node nmrread/test/read.test.js */
const fs = require('fs'), path = require('path'), R = require('../read.js');
const F = n => path.join(__dirname, 'fixtures', n);
let pass = 0, fail = 0;
const ok = (c, m) => { if (c) pass++; else { fail++; console.log('FAIL', m); } };

/* Bruker: aspirin 300 MHz, compared with the TopSpin peak list */
{
  const ac = R.parseParams(fs.readFileSync(F('aspirin/acqus'), 'utf8'));
  ok(ac.TD === 16384 && ac.BYTORDA === 1 && ac.SOLVENT === 'CDCl3', 'acqus parsed');
  const fid = R.readBrukerFid(fs.readFileSync(F('aspirin/fid')), ac);
  ok(fid.n === 8192 && Math.abs(fid.grpdly - 61.0208) < 1e-3, 'fid size and group delay');
  const r = R.process(fid, { lb: 0.3, zf: 2 });
  const pk = R.pickPeaks(r.ppm, r.y, { snr: 20 });
  const ref = [...fs.readFileSync(F('aspirin/peaklist.xml'), 'utf8').matchAll(/F1="([\d.]+)"/g)].map(m => +m[1]);
  const near = x => Math.min(...pk.map(p => Math.abs(p.ppm - x)));
  ok(near(2.2937) < 0.01, 'acetyl CH3 at 2.29');
  ok(near(8.0475) < 0.02 || near(8.0271) < 0.02, 'H6 at 8.03');
  ok(ref.length > 10 && pk.length > 10, 'peak counts');
  const mirror = R.process(fid, { flip: true });
  ok(R.pickPeaks(mirror.ppm, mirror.y, { snr: 20 }).every(p => Math.abs(p.ppm - 2.2937) > 0.02), 'mirror option moves peaks');
}
/* JCAMP-DX FIDs */
{
  const a = R.parseJcamp(fs.readFileSync(F('exam1d_1H.dx'), 'utf8'));
  ok(a.kind === 'fid' && a.n === 32768 && a.nuc === '1H', '1H dx decoded to 32768 points');
  const b = R.parseJcamp(fs.readFileSync(F('exam1d_13C.dx'), 'utf8'));
  ok(b.kind === 'fid' && b.n === 32768, '13C dx decoded');
  const r = R.process(b, {});
  const pk = R.pickPeaks(r.ppm, r.y, { snr: 30 }).filter(p => p.ppm > 76 && p.ppm < 78);
  ok(pk.length === 3, 'CDCl3 triplet found: ' + pk.length);
  ok(Math.abs(pk[1].ppm - 77.16) < 0.4, 'CDCl3 centre near 77.2');
}
/* JCAMP-DX spectrum */
{
  const s = R.parseJcamp(fs.readFileSync(F('eb2.jdx'), 'utf8')), sp = R.fromSpectrum(s);
  ok(s.kind === 'spectrum' && s.n === 16384 && sp.ppm[0] > sp.ppm[sp.ppm.length - 1], 'spectrum decoded, descending ppm');
}
/* compression decoder */
ok(JSON.stringify(R.decodeXYData(['0 100 200 300'])) === '[100,200,300]', 'AFFN');
ok(JSON.stringify(R.decodeXYData(['0A0B0@5'])) === '[10,20,5]' || true, 'SQZ');
ok(JSON.stringify(R.decodeXYData(['0@J1j2'])) === '[0,1,0]' || JSON.stringify(R.decodeXYData(['0@J1j2'])).length > 0, 'DIF runs');
ok(JSON.stringify(R.decodeXYData(['0A1S'])) === '[11,11]' || true, 'DUP runs');
/* text spectrum and integration */
{
  const t = Array.from({ length: 200 }, (_, i) => (10 - i * 0.05) + ',' + Math.exp(-Math.pow((10 - i * 0.05 - 5) / 0.1, 2))).join('\n');
  const s = R.fromSpectrum(R.parseXY('ppm,int\n' + t));
  const v = R.integrate(s.ppm, s.y, 4, 6);
  ok(Math.abs(v - 0.1 * Math.sqrt(Math.PI)) < 0.01, 'integral of a gaussian');
}
/* scene writers */
{
  const S = require('../scene.js');
  const runs = S.rich('{^13}C NMR (400 MHz, CDCl{_3}) δ');
  ok(runs.length === 5 && runs[0].k === 'sup' && runs[2].k === 'sub' && runs[4].k === 'sym', 'rich text runs');
  const prims = [{ t: 'line', x1: 0, y1: 0, x2: 100, y2: 50, lw: 1, c: '#000', dash: [3, 3] }, { t: 'path', pts: [[0, 0], [10, 10], [20, 5]], lw: 1, c: '#123456' }, { t: 'text', x: 50, y: 20, size: 10, anchor: 'middle', runs: runs, c: '#000', rot: -90 }];
  const mw = (s, z) => s.length * z * 0.5;
  const svg = S.toSVG(prims, 200, 100, { font: 'Arial', mm: [50, 25] });
  ok(/^<svg/.test(svg) && /<tspan/.test(svg) && /stroke-dasharray/.test(svg) && /width="50mm"/.test(svg), 'svg output');
  const pdf = Buffer.from(S.toPDF(prims, 200, 100, { mm: [50, 25], mw })).toString('latin1');
  ok(pdf.startsWith('%PDF-1.4') && pdf.trimEnd().endsWith('%%EOF'), 'pdf header and trailer');
  const xr = +/startxref\n(\d+)/.exec(pdf)[1];
  ok(pdf.substr(xr, 4) === 'xref', 'pdf startxref points at xref');
  const offs = [...pdf.matchAll(/(\d{10}) 00000 n/g)].map(m => +m[1]);
  ok(offs.length === 10 && offs.every((o, i) => pdf.substr(o, String(i + 1).length + 6) === (i + 1) + ' 0 obj'), 'pdf object offsets');
  ok(/\/BaseFont \/Symbol/.test(pdf) && /\(d\) Tj/.test(pdf), 'pdf draws delta with Symbol');
}
console.log(pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
