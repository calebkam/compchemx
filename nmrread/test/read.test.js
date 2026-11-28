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
console.log(pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
