/* CCX NMRRead - core reader and processor (pure functions, no DOM).
   Reads Bruker FID directories (acqus + fid), JCAMP-DX (FID or spectrum) and
   plain CSV/text spectra; processes FIDs (digital filter removal, apodisation,
   zero filling, FFT, phasing, baseline, referencing), picks peaks, integrates. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.NMRRead = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  /* ---------- Bruker digital filter table (DSPFVS 10-13, by DECIM) ---------- */
  const DSP = {
    10: { 2: 44.75, 3: 33.5, 4: 66.625, 6: 59.0833333333333, 8: 68.5625, 12: 60.375, 16: 69.53125, 24: 61.0208333333333, 32: 70.015625, 48: 61.34375, 64: 70.2578125, 96: 61.5052083333333, 128: 70.37890625, 192: 61.5859375, 256: 70.439453125, 384: 61.6263020833333, 512: 70.4697265625, 768: 61.6464843750, 1024: 70.48486328125, 1536: 61.6565755208333, 2048: 70.49243164062 },
    11: { 2: 46.0, 3: 36.5, 4: 48.0, 6: 50.1666666666667, 8: 53.25, 12: 69.5, 16: 72.25, 24: 70.1666666666667, 32: 72.75, 48: 70.5, 64: 73.0, 96: 70.6666666666667, 128: 72.5, 192: 71.3333333333333, 256: 72.25, 384: 71.6666666666667, 512: 72.125, 768: 71.8333333333333, 1024: 72.0625, 1536: 71.9166666666667, 2048: 72.03125 },
    12: { 2: 46.0, 3: 36.5, 4: 48.0, 6: 50.1666666666667, 8: 53.25, 12: 69.5, 16: 71.625, 24: 70.1666666666667, 32: 72.125, 48: 70.5, 64: 72.75, 96: 70.6666666666667, 128: 72.5, 192: 71.3333333333333, 256: 72.25, 384: 71.6666666666667, 512: 72.125, 768: 71.8333333333333, 1024: 72.0625, 1536: 71.9166666666667, 2048: 72.03125 },
    13: { 2: 2.75, 3: 2.8333333333333335, 4: 2.875, 6: 2.9166666666666665, 8: 2.9375, 12: 2.9583333333333335, 16: 2.96875, 24: 2.9791666666666665, 32: 2.984375, 48: 2.9895833333333335, 64: 2.9921875, 96: 2.9947916666666665 }
  };

  /* ---------- parameter text (JCAMP-style "##$KEY= value") ---------- */
  function parseParams(text) {
    const out = {}, re = /^##\$?([A-Za-z0-9_.\- ]+?)=\s*([^\r\n]*)/gm;
    let m;
    while ((m = re.exec(text))) {
      let k = m[1].trim(), v = m[2].trim();
      if (k.charAt(0) === '$') k = k.slice(1);
      if (/^<.*>$/.test(v)) out[k] = v.slice(1, -1);
      else if (v !== '' && /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/.test(v)) out[k] = parseFloat(v);
      else out[k] = v;
    }
    return out;
  }

  function groupDelay(p) {
    if (p.GRPDLY > 0) return p.GRPDLY;
    const v = p.DSPFVS, d = p.DECIM;
    if (v === undefined || v >= 14 || v < 10) return 0;
    return (DSP[v] && DSP[v][d]) || 0;
  }

  /* ---------- FFT (radix-2, in place) ---------- */
  function fft(re, im, inverse) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = 2 * Math.PI / len * (inverse ? 1 : -1), wr = Math.cos(ang), wi = Math.sin(ang), h = len >> 1;
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < h; k++) {
          const a = i + k, b = a + h, tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
          re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti;
          const nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
        }
      }
    }
  }
  const pow2 = n => { let p = 1; while (p < n) p <<= 1; return p; };

  /* ---------- Bruker ---------- */
  /* buf: ArrayBuffer/Uint8Array of the fid file; acqus: parsed acqus params */
  function readBrukerFid(buf, acqus) {
    const u8 = buf instanceof Uint8Array ? buf : new Uint8Array(buf);
    const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
    const dbl = acqus.DTYPA === 2, little = acqus.BYTORDA === 0, w = dbl ? 8 : 4;
    const total = Math.floor(u8.byteLength / w), n = Math.floor(total / 2);
    const re = new Float64Array(n), im = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const o = 2 * i * w;
      re[i] = dbl ? dv.getFloat64(o, little) : dv.getInt32(o, little);
      im[i] = dbl ? dv.getFloat64(o + w, little) : dv.getInt32(o + w, little);
    }
    return fidFromBruker(re, im, acqus);
  }
  function fidFromBruker(re, im, p) {
    const n0 = re.length;
    return {
      kind: 'fid', re, im, n: n0,
      sw: p.SW_h, sfo1: p.SFO1 || p.BF1, bf1: p.BF1 || p.SFO1, o1: p.O1 || 0,
      nuc: (p.NUC1 || '1H').replace(/[<>]/g, ''), solvent: p.SOLVENT || '', pulprog: p.PULPROG || '',
      ns: p.NS, aqmod: p.AQ_mod, grpdly: groupDelay(p), params: p, filter: true
    };
  }

  /* remove the Bruker digital filter (nmrglue-compatible, integer group delay) */
  function removeDigitalFilter(fid) {
    if (!fid.filter || !(fid.grpdly > 0)) return fid;
    const n = fid.n, ph = Math.floor(fid.grpdly), skip = Math.floor(ph + 2), add = Math.max(skip - 6, 0);
    if (skip >= n) return fid;
    const re = new Float64Array(n), im = new Float64Array(n);
    for (let m = 0; m < n; m++) { const s = (m + ph) % n; re[m] = fid.re[s]; im[m] = fid.im[s]; }
    for (let i = 0; i < add; i++) { re[i] += re[n - 1 - i]; im[i] += im[n - 1 - i]; }
    const m = n - skip, out = Object.assign({}, fid, { re: re.slice(0, m), im: im.slice(0, m), n: m, filter: false });
    return out;
  }

  /* ---------- JCAMP-DX ---------- */
  const SQZ = '@ABCDEFGHI', SQZn = 'abcdefghi', DIFp = '%JKLMNOPQR', DIFn = 'jklmnopqr', DUPc = 'STUVWXYZs';
  function decodeLine(line, st) {
    /* returns array of y values; st carries {last, diff, lastDif, pendDup} between lines */
    const vals = [];
    const s = line.replace(/\$\$.*$/, '').trim();
    if (!s) return vals;
    let i = 0, first = true, skipCheck = st.skip;
    const L = s.length;
    /* X value (first token, up to first SQZ/DIF letter or sign after digits) */
    const m0 = /^[+-]?(\d+\.?\d*|\.\d+)/.exec(s);
    if (m0) { i = m0[0].length; st.x = parseFloat(m0[0]); }
    st.lastDif = false;
    let lastDifToken = false;
    while (i < L) {
      const c = s[i];
      if (c === ' ' || c === ',' || c === '\t') { i++; continue; }
      let mode = 'a', sign = 1, lead = '', dupN = 0;
      if (SQZ.indexOf(c) >= 0) { mode = 's'; lead = String(SQZ.indexOf(c)); i++; }
      else if (SQZn.indexOf(c) >= 0) { mode = 's'; sign = -1; lead = String(SQZn.indexOf(c) + 1); i++; }
      else if (DIFp.indexOf(c) >= 0) { mode = 'd'; lead = String(DIFp.indexOf(c)); i++; }
      else if (DIFn.indexOf(c) >= 0) { mode = 'd'; sign = -1; lead = String(DIFn.indexOf(c) + 1); i++; }
      else if (DUPc.indexOf(c) >= 0) { mode = 'u'; lead = String(DUPc.indexOf(c) + 1); i++; }
      else if (c === '+' || c === '-') { mode = 'a'; if (c === '-') sign = -1; i++; }
      else if (/[0-9.]/.test(c)) { mode = 'a'; }
      else { i++; continue; }
      let j = i;
      if (mode === 'a') { while (j < L && /[0-9.eE]/.test(s[j]) && !(s[j] === 'e' || s[j] === 'E' ? !/[+-]?\d/.test(s.slice(j + 1, j + 3)) : false)) { if ((s[j] === 'e' || s[j] === 'E')) { j++; if (s[j] === '+' || s[j] === '-') j++; } else j++; } }
      else { while (j < L && /[0-9]/.test(s[j])) j++; }
      const digits = lead + s.slice(i, j);
      i = j;
      if (mode === 'u') {
        const cnt = parseInt(digits, 10);
        for (let r = 1; r < cnt; r++) pushVal(st.lastTok);
        continue;
      }
      const num = sign * parseFloat(digits || '0');
      st.lastTok = { mode, num };
      pushVal(st.lastTok);
      lastDifToken = mode === 'd';
    }
    function pushVal(t) {
      if (t.mode === 'd') { st.last += t.num; st.diff = true; if (skipCheck && first) { /* never happens: DIF can't start a line */ } vals.push(st.last); }
      else {
        if (skipCheck && first) { first = false; skipCheck = false; st.last = t.num; return; } /* DIF check value */
        first = false; st.last = t.num; st.diff = false; vals.push(t.num);
      }
    }
    st.lastDif = lastDifToken;
    return vals;
  }
  function decodeXYData(lines) {
    const st = { last: 0, diff: false, lastDif: false, lastTok: null }, out = [];
    for (const ln of lines) {
      /* a DIF line's last value is repeated as the first value of the next line; the X index tells which case applies */
      const mx = /^\s*([+-]?\d+)(?![\d.eE])/.exec(ln.replace(/\$\$.*$/, ''));
      st.skip = st.lastDif;
      if (mx && st.lastDif) { const x = parseInt(mx[1], 10); if (x === out.length - 1) st.skip = true; else if (x === out.length) st.skip = false; }
      const v = decodeLine(ln, st); for (let k = 0; k < v.length; k++) out.push(v[k]);
    }
    return out;
  }

  function splitBlocks(text) {
    const lines = text.split(/\r?\n/), blocks = [];
    let cur = null;
    for (const ln of lines) {
      if (/^##TITLE=/i.test(ln)) { cur = { head: [], lines: [] }; blocks.push(cur); }
      if (!cur) { cur = { head: [], lines: [] }; blocks.push(cur); }
      cur.lines.push(ln);
    }
    return blocks;
  }

  function parseJcamp(text) {
    const params = parseParams(text), blocks = splitBlocks(text), found = [];
    const meta = k => { const m = new RegExp('^##' + k + '=\\s*(.*)$', 'im').exec(text); return m ? m[1].trim() : ''; };
    for (const b of blocks) {
      const lines = b.lines, ntIdx = lines.findIndex(l => /^##NTUPLES=/i.test(l));
      if (ntIdx >= 0) {
        const nt = (lines[ntIdx].split('=')[1] || '').trim().toUpperCase();
        const arr = k => { const l = lines.find(x => new RegExp('^##' + k + '=', 'i').test(x)); return l ? l.split('=').slice(1).join('=').split(',').map(s => s.trim()) : []; };
        const names = arr('VAR_NAME'), fac = arr('FACTOR').map(Number), first = arr('FIRST').map(Number), last = arr('LAST').map(Number), dim = arr('VAR_DIM').map(Number), units = arr('UNITS');
        const pages = [];
        let page = null;
        for (let i = ntIdx; i < lines.length; i++) {
          const l = lines[i];
          if (/^##PAGE=/i.test(l)) { page = { n: l, data: [] }; pages.push(page); continue; }
          if (/^##DATA TABLE=/i.test(l)) continue;
          if (/^##END/i.test(l)) { page = null; continue; }
          if (/^##/.test(l)) { if (page) page = null; continue; }
          if (page && l.trim() && !/^\$\$/.test(l)) page.data.push(l);
        }
        const dec = pages.map(p => decodeXYData(p.data));
        const isFid = /FID/.test(nt) || names.some(n => /FID/i.test(n));
        const re = dec[0], im = dec[1] || null;
        if (!re || !re.length) continue;
        const fR = fac[1] || 1, fI = fac[2] || 1;
        const R = Float64Array.from(re, v => v * fR), I = im ? Float64Array.from(im, v => v * fI) : null;
        found.push({ isFid, re: R, im: I, xFirst: first[0], xLast: last[0], xUnits: (units[0] || '').toUpperCase(), n: R.length });
        continue;
      }
      const xyIdx = lines.findIndex(l => /^##(XYDATA|PEAK TABLE|XYPOINTS)=/i.test(l));
      if (xyIdx >= 0) {
        const get = k => { const l = lines.find(x => new RegExp('^##' + k + '=', 'i').test(x)); return l ? parseFloat(l.split('=')[1]) : NaN; };
        const data = [];
        for (let i = xyIdx + 1; i < lines.length; i++) { if (/^##/.test(lines[i])) break; if (lines[i].trim()) data.push(lines[i]); }
        const y = decodeXYData(data), fy = get('YFACTOR') || 1, fx = get('XFACTOR') || 1;
        found.push({ isFid: false, re: Float64Array.from(y, v => v * fy), im: null, xFirst: get('FIRSTX') * fx, xLast: get('LASTX') * fx, xUnits: (meta('XUNITS') || '').toUpperCase(), n: y.length, simple: true });
      }
    }
    const pick = found.find(f => f.isFid) || found[0];
    if (!pick) throw new Error('No data found in this JCAMP-DX file.');
    const info = {
      title: meta('TITLE'), solvent: params.SOLVENT || meta('\\.SOLVENT NAME') || '', nuc: (params.NUC1 || meta('\\.OBSERVE NUCLEUS') || '1H').replace(/[<>^]/g, ''),
      sfo1: params.SFO1 || parseFloat(meta('\\.OBSERVE FREQUENCY')) || params.BF1, bf1: params.BF1 || params.SFO1 || parseFloat(meta('\\.OBSERVE FREQUENCY')), o1: params.O1 || 0
    };
    if (pick.isFid) {
      const sw = params.SW_h || (pick.xLast && pick.n > 1 ? pick.n / (pick.xLast + (pick.xLast / (pick.n - 1))) : 0);
      const fid = {
        kind: 'fid', re: pick.re, im: pick.im || new Float64Array(pick.n), n: pick.n, sw: sw, sfo1: info.sfo1, bf1: info.bf1, o1: info.o1, nuc: info.nuc,
        solvent: info.solvent, pulprog: params.PULPROG || '', ns: params.NS, aqmod: params.AQ_mod, grpdly: groupDelay(params), params, title: info.title, filter: pick.n * 2 === params.TD ? true : true
      };
      /* JCAMP FIDs exported from TopSpin hold the raw (filtered) points */
      return fid;
    }
    /* processed spectrum: build ppm axis */
    const sf = info.sfo1 || info.bf1 || 400, xs = pick.xUnits === 'PPM' ? [pick.xFirst, pick.xLast] : [pick.xFirst / sf, pick.xLast / sf];
    return { kind: 'spectrum', re: pick.re, im: pick.im, n: pick.n, ppmFirst: xs[0], ppmLast: xs[1], sfo1: sf, nuc: info.nuc, solvent: info.solvent, title: info.title, params };
  }

  /* ---------- plain text / CSV spectra: "ppm, intensity" ---------- */
  function parseXY(text) {
    const xs = [], ys = [];
    for (const l of text.split(/\r?\n/)) {
      if (/^\s*(#|$|[A-Za-z])/.test(l)) continue;
      const t = l.trim().split(/[\s,;\t]+/).map(Number);
      if (t.length >= 2 && isFinite(t[0]) && isFinite(t[1])) { xs.push(t[0]); ys.push(t[1]); }
    }
    if (xs.length < 8) throw new Error('Could not read two numeric columns (ppm, intensity).');
    return { kind: 'spectrum', re: Float64Array.from(ys), im: null, n: ys.length, ppmFirst: xs[0], ppmLast: xs[xs.length - 1], sfo1: 400, nuc: '1H', solvent: '', title: '', params: {}, xs: Float64Array.from(xs) };
  }

  /* ---------- processing ---------- */
  const DEF = { lb: 0.3, zf: 2, p0: 0, p1: 0, pivot: 0.5, auto: true, flip: false, refShift: 0, base: true, baseOrder: 3, rev: false };

  /* time-domain: apodise + zero-fill + FFT; returns complex spectrum in ppm order high->low */
  function transform(fid0, o) {
    o = Object.assign({}, DEF, o);
    const fid = o.keepFilter ? fid0 : removeDigitalFilter(fid0);
    const n = fid.n, size = pow2(Math.max(n, Math.round(n * o.zf)));
    const re = new Float64Array(size), im = new Float64Array(size), dt = 1 / fid.sw;
    for (let i = 0; i < n; i++) {
      const w = Math.exp(-Math.PI * o.lb * i * dt) * (o.gm ? Math.exp(-Math.pow(i * dt / (o.gm), 2) * 0) : 1);
      re[i] = fid.re[i] * w * (i === 0 ? 0.5 : 1); im[i] = fid.im[i] * w * (i === 0 ? 0.5 : 1);
    }
    fft(re, im, false);
    /* fftshift, then orient so ppm decreases left->right */
    const half = size >> 1, R = new Float64Array(size), I = new Float64Array(size);
    for (let k = 0; k < size; k++) { const s = (k + half) % size; R[k] = re[s]; I[k] = im[s]; }
    if (!o.flip) { R.reverse(); I.reverse(); }
    /* axis: centre = carrier, in ppm */
    const sf = fid.sfo1 || fid.bf1 || 400, bf = fid.bf1 || sf, swp = fid.sw / sf, centre = (fid.o1 || 0) / bf + (o.refShift || 0);
    const ppm = new Float64Array(size);
    for (let k = 0; k < size; k++) ppm[k] = centre + swp / 2 - swp * k / size;
    return { re: R, im: I, ppm, sf, sw: fid.sw, size, fid };
  }

  function phased(S, p0, p1, pivot) {
    const n = S.re.length, re = new Float64Array(n), im = new Float64Array(n), a0 = p0 * Math.PI / 180, a1 = p1 * Math.PI / 180, pv = (pivot === undefined ? 0.5 : pivot);
    for (let k = 0; k < n; k++) {
      const a = a0 + a1 * (k / n - pv), c = Math.cos(a), s = Math.sin(a);
      re[k] = S.re[k] * c - S.im[k] * s; im[k] = S.re[k] * s + S.im[k] * c;
    }
    return { re, im };
  }

  /* entropy-minimisation auto phase (ACME) */
  function autoPhase(S) {
    const n = S.re.length;
    const cost = (p0, p1) => {
      const P = phased(S, p0, p1, 0.5), re = P.re;
      let sumAbs = 0, d = new Float64Array(n - 1);
      for (let k = 0; k < n - 1; k++) { d[k] = Math.abs(re[k + 1] - re[k]); sumAbs += d[k]; }
      if (!sumAbs) return 0;
      let h = 0, pen = 0;
      for (let k = 0; k < n - 1; k++) { const p = d[k] / sumAbs; if (p > 0) h -= p * Math.log(p); }
      let negSq = 0;
      for (let k = 0; k < n; k++) if (re[k] < 0) negSq += re[k] * re[k];
      let tot = 0; for (let k = 0; k < n; k++) tot += re[k] * re[k];
      pen = tot ? 1000 * negSq / tot : 0;
      return h + pen;
    };
    let best = [0, 0], bc = Infinity;
    for (let a = -180; a < 180; a += 6) { const c = cost(a, 0); if (c < bc) { bc = c; best = [a, 0]; } }
    let step0 = 6, step1 = 40;
    for (let it = 0; it < 24; it++) {
      let moved = false;
      for (const [d0, d1] of [[step0, 0], [-step0, 0], [0, step1], [0, -step1]]) {
        const c = cost(best[0] + d0, best[1] + d1);
        if (c < bc - 1e-12) { bc = c; best = [best[0] + d0, best[1] + d1]; moved = true; }
      }
      if (!moved) { step0 /= 2; step1 /= 2; if (step0 < 0.05) break; }
    }
    let p0 = ((best[0] + 180) % 360 + 360) % 360 - 180;
    return { p0, p1: best[1] };
  }

  /* baseline: iterative polynomial fit to points below a noise-based threshold */
  function baselineCorrect(y, order) {
    const n = y.length, m = Math.min(order, 8);
    const x = new Float64Array(n); for (let i = 0; i < n; i++) x[i] = 2 * i / (n - 1) - 1;
    let mask = new Uint8Array(n).fill(1), coef = null;
    const fit = () => {
      const A = Array.from({ length: m + 1 }, () => new Float64Array(m + 2));
      for (let i = 0; i < n; i++) if (mask[i]) {
        const pw = [1]; for (let j = 1; j <= 2 * m; j++) pw[j] = pw[j - 1] * x[i];
        for (let r = 0; r <= m; r++) { for (let c = 0; c <= m; c++) A[r][c] += pw[r + c]; A[r][m + 1] += pw[r] * y[i]; }
      }
      for (let c = 0; c <= m; c++) {
        let p = c; for (let r = c + 1; r <= m; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
        [A[c], A[p]] = [A[p], A[c]];
        if (!A[c][c]) continue;
        for (let r = 0; r <= m; r++) if (r !== c) { const f = A[r][c] / A[c][c]; for (let k = c; k <= m + 1; k++) A[r][k] -= f * A[c][k]; }
      }
      return A.map((row, i) => row[m + 1] / row[i]);
    };
    const ev = (c, xv) => { let s = 0, p = 1; for (let j = 0; j <= m; j++) { s += c[j] * p; p *= xv; } return s; };
    for (let it = 0; it < 12; it++) {
      coef = fit();
      let ss = 0, cnt = 0; const res = new Float64Array(n);
      for (let i = 0; i < n; i++) { res[i] = y[i] - ev(coef, x[i]); if (mask[i]) { ss += res[i] * res[i]; cnt++; } }
      const sd = Math.sqrt(ss / Math.max(cnt, 1));
      for (let i = 0; i < n; i++) mask[i] = res[i] < 2.5 * sd && res[i] > -4 * sd ? 1 : 0;
    }
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) out[i] = y[i] - ev(coef, x[i]);
    return out;
  }

  /* full pipeline from a fid to a real spectrum */
  function process(fid, o) {
    o = Object.assign({}, DEF, o);
    const S = transform(fid, o);
    let p0 = o.p0, p1 = o.p1;
    if (o.auto) { const a = autoPhase(S); p0 = a.p0; p1 = a.p1; }
    const P = phased(S, p0, p1, 0.5);
    let y = P.re;
    /* make the tallest feature positive */
    let mx = 0, mn = 0; for (let k = 0; k < y.length; k++) { if (y[k] > mx) mx = y[k]; if (y[k] < mn) mn = y[k]; }
    let flipSign = false;
    if (o.auto && -mn > mx) { y = y.map(v => -v); flipSign = true; p0 += 180; }
    if (o.base) y = baselineCorrect(y, o.baseOrder);
    return { ppm: S.ppm, y: Float64Array.from(y), p0: ((p0 + 180) % 360 + 360) % 360 - 180, p1, sf: S.sf, sw: S.sw, S, flipSign };
  }

  /* from a processed-spectrum import */
  function fromSpectrum(sp) {
    const n = sp.n, ppm = new Float64Array(n);
    if (sp.xs) for (let i = 0; i < n; i++) ppm[i] = sp.xs[i];
    else for (let i = 0; i < n; i++) ppm[i] = sp.ppmFirst + (sp.ppmLast - sp.ppmFirst) * i / (n - 1);
    let y = Float64Array.from(sp.re);
    if (ppm[0] < ppm[n - 1]) { ppm.reverse(); y.reverse(); }
    return { ppm, y, sf: sp.sfo1 || 400, p0: 0, p1: 0 };
  }

  /* ---------- noise, peaks, integrals ---------- */
  function noiseLevel(y) {
    const n = y.length, seg = Math.max(64, n >> 6), sds = [];
    for (let s = 0; s + seg <= n; s += seg) {
      let m = 0; for (let i = s; i < s + seg; i++) m += y[i]; m /= seg;
      let v = 0; for (let i = s; i < s + seg; i++) v += (y[i] - m) * (y[i] - m);
      sds.push(Math.sqrt(v / seg));
    }
    sds.sort((a, b) => a - b);
    return sds[Math.floor(sds.length * 0.1)] || sds[0] || 1;
  }
  function pickPeaks(ppm, y, o) {
    o = Object.assign({ snr: 8, lo: -Infinity, hi: Infinity }, o);
    const thr = (o.abs || 0) || o.snr * noiseLevel(y), out = [];
    for (let i = 1; i < y.length - 1; i++) {
      if (y[i] > thr && y[i] >= y[i - 1] && y[i] > y[i + 1] && ppm[i] <= o.hi && ppm[i] >= o.lo) {
        const a = y[i - 1], b = y[i], c = y[i + 1], den = a - 2 * b + c, d = den ? 0.5 * (a - c) / den : 0;
        out.push({ ppm: ppm[i] + d * (ppm[i + 1] - ppm[i]), h: b - 0.25 * (a - c) * d, i });
      }
    }
    return out;
  }
  function integrate(ppm, y, lo, hi) {
    let s = 0;
    for (let i = 1; i < y.length; i++) if (ppm[i] <= hi && ppm[i] >= lo) s += y[i] * Math.abs(ppm[i] - ppm[i - 1]);
    return s;
  }

  /* residual solvent shifts for quick referencing (1H, ppm) */
  const SOLVENT_1H = { 'CDCl3': 7.26, 'DMSO': 2.50, 'DMSO-D6': 2.50, 'D2O': 4.79, 'CD3OD': 3.31, 'ACETONE': 2.05, 'C6D6': 7.16, 'CD3CN': 1.94, 'CD2Cl2': 5.32, 'TMS': 0 };
  const SOLVENT_13C = { 'CDCl3': 77.16, 'DMSO': 39.52, 'DMSO-D6': 39.52, 'CD3OD': 49.00, 'ACETONE': 29.84, 'C6D6': 128.06, 'CD3CN': 1.32, 'CD2Cl2': 53.84 };

  /* ---------- file entry point ---------- */
  function readAny(name, content) {
    if (/\.(jdx|dx|jcm|jcamp)$/i.test(name) || (typeof content === 'string' && /^##TITLE=/m.test(content.slice(0, 4000)))) return parseJcamp(typeof content === 'string' ? content : new TextDecoder().decode(content));
    if (/\.(csv|txt|tsv|xy)$/i.test(name)) return parseXY(typeof content === 'string' ? content : new TextDecoder().decode(content));
    throw new Error('Unsupported file type: ' + name);
  }

  return { DSP, parseParams, groupDelay, fft, readBrukerFid, fidFromBruker, removeDigitalFilter, parseJcamp, parseXY, decodeXYData, readAny, transform, phased, autoPhase, baselineCorrect, process, fromSpectrum, noiseLevel, pickPeaks, integrate, SOLVENT_1H, SOLVENT_13C, DEF };
});
