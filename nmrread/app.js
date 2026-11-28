/* CCX NMRRead - page logic: load, process, view, integrate, export. */
(function () {
  'use strict';
  const R = window.NMRRead, $ = id => document.getElementById(id);
  const S = { fid: null, ppm0: null, y: null, ppm: null, view: null, vz: 1, ints: [], scale: null, refShift: 0, peaks: [], rawPeaks: [], exps: [], masks: [], mode: null, yRaw: null, refMode: false, nuc: '1H', meta: {}, kind: '' };
  const stage = $('stage'); let svgEl = null, lay = null, raf = 0;

  /* ---------- helpers ---------- */
  const SUP = { 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' }, SUB = { 0: '₀', 1: '₁', 2: '₂', 3: '₃', 4: '₄', 5: '₅', 6: '₆', 7: '₇', 8: '₈', 9: '₉' };
  const nucLabel = n => (n || '1H').replace(/^\d+/, d => d.replace(/\d/g, c => SUP[c]));
  const solvLabel = s => (s || '').replace(/(\D)(\d)/g, (m, a, d) => a + SUB[d]).replace(/-d(\d)/i, '-d$1');
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  function toast(t) { const e = $('toast'); e.textContent = t; e.hidden = false; clearTimeout(toast.t); toast.t = setTimeout(() => e.hidden = true, 3200); }
  const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const mhz = () => { const f = S.fid ? (S.fid.bf1 || S.fid.sfo1) : (S.meta.sf || 0); return f ? Math.round(f) : null; };

  /* ---------- loading ---------- */
  async function walk(entry, path, out, depth) {
    if (out.length > 600 || depth > 6) return;
    if (entry.isFile) { await new Promise(res => entry.file(f => { out.push({ path: path + entry.name, file: f }); res(); }, res)); return; }
    if (!entry.isDirectory) return;
    const rd = entry.createReader();
    for (;;) {
      const batch = await new Promise(res => rd.readEntries(res, () => res([])));
      if (!batch.length) break;
      for (const e of batch) await walk(e, path + entry.name + '/', out, depth + 1);
    }
  }
  async function fromDrop(dt) {
    const out = [];
    if (dt.items && dt.items.length && dt.items[0].webkitGetAsEntry) {
      const ents = []; for (const it of dt.items) { const e = it.webkitGetAsEntry(); if (e) ents.push(e); }
      for (const e of ents) await walk(e, '', out, 0);
    } else for (const f of dt.files) out.push({ path: f.name, file: f });
    return out;
  }
  const dirOf = p => p.replace(/[^/]*$/, '');
  async function load(list) {
    try {
      const lc = list.map(x => Object.assign(x, { n: x.file.name.toLowerCase() }));
      if (lc.some(x => x.n === 'ser') && !lc.some(x => x.n === 'fid')) throw new Error('This looks like 2D data (ser). NMRRead reads 1D data only.');
      const fids = lc.filter(x => x.n === 'fid').sort((a, b) => a.path.length - b.path.length);
      let fid = null, title = '';
      if (fids.length) {
        const f = fids[0], dir = dirOf(f.path), aq = lc.find(x => x.n === 'acqus' && dirOf(x.path) === dir) || lc.find(x => x.n === 'acqus');
        if (!aq) throw new Error('Found a fid but no acqus file next to it. Open the whole experiment folder, or add acqus.');
        const p = R.parseParams(await aq.file.text());
        if (p.SW_h === undefined || !p.TD) throw new Error('acqus is missing SW_h or TD.');
        fid = R.readBrukerFid(await f.file.arrayBuffer(), p);
        title = f.path;
        if (fids.length > 1) toast('Opened ' + f.path + ' (' + fids.length + ' experiments found)');
        S.meta = {};
      } else {
        const f = lc.find(x => /\.(jdx|dx|jcm|jcamp|csv|txt|tsv|xy)$/.test(x.n));
        if (!f) throw new Error('No readable data found. Use a Bruker folder (fid + acqus), a .jdx/.dx file, or a two-column text file.');
        const d = R.readAny(f.file.name, await f.file.text());
        title = f.file.name;
        if (d.kind === 'fid') fid = d; else { setSpectrum(d, title); return; }
      }
      S.fid = fid; S.kind = 'fid'; S.title = title; S.refShift = 0; S.ints = []; S.exps = []; S.masks = []; S.scale = null; S.vz = 1;
      $('b-auto').classList.add('on'); S.auto = true; $('flip').checked = false;
      S.nuc = fid.nuc;
      fillSolvents(); preselectSolvent(fid.solvent);
      setControls(true);
      recompute(true);
    } catch (e) { toast(e.message || String(e)); console.error(e); }
  }
  function setSpectrum(d, title) {
    const sp = R.fromSpectrum(d);
    S.fid = null; S.kind = 'spectrum'; S.title = title; S.ppm0 = sp.ppm; S.yRaw = sp.y; S.nuc = d.nuc; S.refShift = 0; S.ints = []; S.exps = []; S.masks = []; S.scale = null; S.vz = 1;
    S.meta = { sf: d.sfo1, solvent: d.solvent, nuc: d.nuc };
    fillSolvents(); preselectSolvent(d.solvent);
    setControls(false);
    finishSpectrum(true);
  }
  function setControls(isFid) {
    for (const id of ['lb', 'zf', 'b-auto', 'p0', 'p1', 'base', 'bord', 'flip']) $(id).disabled = !isFid;
  }

  /* ---------- processing ---------- */
  function recompute(resetView) {
    if (!S.fid) return;
    const o = { lb: +$('lb').value, zf: +$('zf').value, auto: S.auto, p0: +$('p0').value, p1: +$('p1').value, flip: $('flip').checked, base: $('base').checked, baseOrder: +$('bord').value };
    const r = R.process(S.fid, o);
    S.ppm0 = r.ppm; S.yRaw = r.y;
    if (S.auto) { $('p0').value = r.p0; $('p1').value = r.p1; }
    $('p0v').textContent = Math.round(+$('p0').value * 10) / 10 + '°'; $('p1v').textContent = Math.round(+$('p1').value) + '°';
    $('phv').textContent = S.auto ? 'auto' : 'manual';
    $('lbv').textContent = (+$('lb').value).toFixed(2) + ' Hz';
    finishSpectrum(resetView);
  }
  function finishSpectrum(resetView) {
    const n = S.ppm0.length;
    if (!S.ppm || S.ppm.length !== n) S.ppm = new Float64Array(n);
    for (let i = 0; i < n; i++) S.ppm[i] = S.ppm0[i] + S.refShift;
    S.rawPeaks = R.pickPeaks(S.ppm, S.yRaw, { snr: +$('snr').value });
    applyMasks();
    S.noise = R.noiseLevel(S.y);
    S.peaks = R.pickPeaks(S.ppm, S.y, { snr: +$('snr').value });
    if (resetView || !S.view) S.view = fullView();
    $('msg').hidden = true;
    renderInfo(); renderTable(); renderChips(); draw();
  }
  function fullView() {
    const p = S.ppm, y = S.y, th = S.noise * 6;
    let hi = -Infinity, lo = Infinity;
    for (let i = 0; i < y.length; i++) if (Math.abs(y[i]) > th) { if (p[i] > hi) hi = p[i]; if (p[i] < lo) lo = p[i]; }
    if (!isFinite(hi) || hi - lo < 0.05) { hi = p[0]; lo = p[p.length - 1]; }
    const pad = Math.max(0.3, (hi - lo) * 0.05) * (S.nuc === '1H' ? 1 : 4);
    return [Math.min(p[0], hi + pad), Math.max(p[p.length - 1], lo - pad)];
  }
  function renderInfo() {
    const f = S.fid, m = S.meta, L = [];
    L.push(S.title || '');
    L.push('Nucleus ' + S.nuc + (mhz() ? ' · ' + mhz() + ' MHz' : ''));
    if (f) { L.push('Points ' + f.n + ' · SW ' + (f.sw / (f.sfo1 || 1)).toFixed(2) + ' ppm'); if (f.grpdly) L.push('Digital filter removed (' + f.grpdly.toFixed(1) + ' points)'); if (f.pulprog) L.push('Pulse programme ' + f.pulprog + (f.ns ? ' · NS ' + f.ns : '')); if (f.solvent) L.push('Solvent ' + f.solvent); }
    else if (m.solvent) L.push('Solvent ' + m.solvent);
    $('info').textContent = L.filter(Boolean).join('\n');
  }

  /* ---------- solvent referencing ---------- */
  function fillSolvents() {
    const t = S.nuc === '13C' ? R.SOLVENT_13C : R.SOLVENT_1H, s = $('solv');
    s.innerHTML = '<option value="">—</option>' + Object.keys(t).filter(k => k !== 'DMSO-D6').map(k => '<option value="' + k + '">' + esc(solvLabel(k === 'DMSO' ? 'DMSO-d6' : k)) + ' (' + t[k] + ')</option>').join('');
  }
  function preselectSolvent(name) {
    if (!name) return;
    const u = name.toUpperCase().replace(/[-_]D\d?$/, '').replace(/-D6$/, '');
    for (const o of $('solv').options) if (o.value && (o.value.toUpperCase() === name.toUpperCase() || o.value.toUpperCase() === u)) { $('solv').value = o.value; return; }
  }
  function setSolvent() {
    const k = $('solv').value; if (!k || !S.ppm) return toast('Choose the solvent first.');
    const t = (S.nuc === '13C' ? R.SOLVENT_13C : R.SOLVENT_1H)[k]; if (t === undefined) return;
    const win = S.nuc === '13C' ? 1.0 : 0.35, c = S.rawPeaks.filter(p => Math.abs(p.ppm - t) < win);
    if (!c.length) return toast('No peak found near the solvent signal. Lower the peak threshold or use Pick peak.');
    const mx = Math.max.apply(null, c.map(p => p.h)), use = c.filter(p => p.h > mx * 0.3);
    const cen = use.reduce((a, p) => a + p.ppm * p.h, 0) / use.reduce((a, p) => a + p.h, 0);
    S.refShift += t - cen; finishSpectrum(false); toast('Referenced to ' + k + ' at ' + t + ' ppm');
    $('b-refclr').hidden = false;
  }

  /* ---------- suppression (solvent, common impurities, manual) ---------- */
  const SOLV_WIN = {
    '1H': { CDCl3: [7.26, .05], DMSO: [2.50, .07], CD3OD: [3.31, .07], D2O: [4.79, .25], ACETONE: [2.05, .06], C6D6: [7.16, .06], CD3CN: [1.94, .06], CD2Cl2: [5.32, .05] },
    '13C': { CDCl3: [77.16, .6], DMSO: [39.52, 1.0], CD3OD: [49.0, 1.0], ACETONE: [29.84, .8], C6D6: [128.06, .6], CD3CN: [1.32, .8], CD2Cl2: [53.84, .6] }
  };
  const WATER = { CDCl3: 1.56, DMSO: 3.33, CD3OD: 4.87, ACETONE: 2.84, C6D6: 0.40, CD3CN: 2.13, CD2Cl2: 1.52 };
  function autoMasks() {
    const out = [], k = $('solv').value, tab = SOLV_WIN[S.nuc === '13C' ? '13C' : '1H'];
    if ($('m-solv').checked && k && tab[k]) out.push({ lo: tab[k][0] - tab[k][1], hi: tab[k][0] + tab[k][1], auto: 1 });
    if ($('m-imp').checked && S.nuc === '1H') {
      if (WATER[k]) out.push({ lo: WATER[k] - 0.12, hi: WATER[k] + 0.12, auto: 1 });
      out.push({ lo: 0.01, hi: 0.13, auto: 1 });
      if (k === 'CDCl3' || k === 'C6D6' || !k) { out.push({ lo: 1.21, hi: 1.31, auto: 1 }); out.push({ lo: 0.83, hi: 0.93, auto: 1 }); }
    }
    return out;
  }
  function applyMasks() {
    const y = Float64Array.from(S.yRaw), p = S.ppm, all = autoMasks().concat(S.masks);
    for (const m of all) {
      const a = bsearch(p, m.hi), b = Math.min(p.length - 1, bsearch(p, m.lo));
      if (b <= a) continue;
      const ya = y[Math.max(0, a - 1)], yb = y[Math.min(p.length - 1, b + 1)];
      for (let i = a; i <= b; i++) y[i] = ya + (yb - ya) * (i - a) / (b - a);
    }
    S.y = y;
  }
  function bsearch(p, val) { /* first index with p[i] <= val (p descending) */
    let a = 0, b = p.length - 1;
    while (a < b) { const m = (a + b) >> 1; if (p[m] <= val) b = m; else a = m + 1; }
    return a;
  }
  const fmt = (a, b) => a.toFixed(2) + '–' + b.toFixed(2);
  function renderChips() {
    const el = $('chips');
    el.innerHTML = S.exps.map((r, i) => '<span class="pill">Expansion ' + fmt(r.hi, r.lo) + '<button class="x" data-e="' + i + '" title="Remove">✕</button></span>').join('') +
      S.masks.map((r, i) => '<span class="pill">Masked ' + fmt(r.hi, r.lo) + '<button class="x" data-m="' + i + '" title="Remove">✕</button></span>').join('');
  }

  /* ---------- scene ---------- */
  const SC = window.NMRScene;
  let measureCtx = null;
  function mw(str, size, bold, font) {
    if (!measureCtx) measureCtx = document.createElement('canvas').getContext('2d');
    measureCtx.font = (bold ? '600 ' : '') + size + 'px ' + (font || mw.font || 'Arial, Helvetica, sans-serif');
    return measureCtx.measureText(str).width;
  }
  function niceStep(span, target) {
    const raw = span / target, pw = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / pw;
    const f = m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10;
    return { step: f * pw, minor: f === 2 ? 4 : 5 };
  }
  function rangeMax(p, y, hi, lo) {
    let m = 0;
    for (let i = Math.max(0, bsearch(p, hi) - 1); i < p.length && p[i] >= lo; i++) if (y[i] > m) m = y[i];
    return m || 1;
  }
  function intVal(r) { return R.integrate(S.ppm, S.y, r.lo, r.hi); }
  function intH(r) { return S.scale ? intVal(r) / S.scale : 0; }
  const nucRich = n => (n || '1H').replace(/^(\d+)/, '{^$1}');
  const solvRich = s => (s || '').replace(/(\D)(\d+)/g, '$1{_$2}');
  function headLine() {
    const m = mhz(), sol = ((S.fid ? S.fid.solvent : S.meta.solvent) || '').replace(/^DMSO$/i, 'DMSO-d6');
    return nucRich(S.nuc) + ' NMR' + (m || sol ? ' (' + (m ? m + ' MHz' : '') + (m && sol ? ', ' : '') + solvRich(sol) + ')' : '');
  }

  /* trace of p/y between hi and lo, mapped by X and a y function, clamped at yTop */
  function tracePts(hi, lo, X, Yf, yTop, cols, xl) {
    const p = S.ppm, y = S.y, n = p.length, i0 = Math.max(0, bsearch(p, hi) - 2), i1 = Math.min(n - 1, bsearch(p, lo) + 2), pts = [];
    const cl = v => Math.max(yTop, v);
    if (i1 - i0 < cols * 3) { for (let i = i0; i <= i1; i++) pts.push([X(p[i]), cl(Yf(y[i]))]); }
    else {
      const per = (i1 - i0 + 1) / cols;
      for (let c = 0; c < cols; c++) {
        const a = Math.floor(i0 + c * per), b = Math.min(i1, Math.floor(i0 + (c + 1) * per));
        let mn = Infinity, mx = -Infinity; for (let i = a; i <= b; i++) { if (y[i] < mn) mn = y[i]; if (y[i] > mx) mx = y[i]; }
        const x = xl + c + 0.5, up = c % 2 === 0;
        pts.push([x, cl(Yf(up ? mn : mx))], [x, cl(Yf(up ? mx : mn))]);
      }
    }
    return pts;
  }
  function axisTicks(prims, X, hi, lo, axisY, fs, lw, fg, labelEvery) {
    const nt = niceStep(hi - lo, labelEvery), minor = nt.step / nt.minor, dec = Math.max(0, -Math.floor(Math.log10(nt.step) + 1e-9));
    for (let t = Math.ceil(lo / minor - 1e-9) * minor; t <= hi + 1e-9; t += minor) {
      const major = Math.abs(t / nt.step - Math.round(t / nt.step)) < 1e-6, x = X(t);
      prims.push({ t: 'line', x1: x, y1: axisY, x2: x, y2: axisY + (major ? fs * 0.55 : fs * 0.3), lw, c: fg });
      if (major) prims.push({ t: 'text', x, y: axisY + fs * 1.65, size: fs, anchor: 'middle', runs: [{ s: (Math.abs(t) < 1e-9 ? 0 : t).toFixed(dec), k: 'n' }], c: fg });
    }
  }

  /* Build the spectrum as drawing primitives. o: {fs, lw, col, fg, labels, ints, curve, info, title} */
  function build(W, H, o) {
    const v = S.view, hi = v[0], lo = v[1], fs = o.fs, lw = o.lw, prims = [];
    const showInt = o.ints && S.ints.length, tickH = fs * 1.9;
    const xl = fs * 1.2, xr = W - fs * 1.2, span = hi - lo, X = p => xl + (hi - p) / span * (xr - xl);
    /* integral labels: first row where they do not touch the previous one */
    const place = []; let nrows = 1;
    if (showInt) {
      const ends = [-1e9, -1e9, -1e9];
      for (const r of S.ints) {
        if (r.hi < lo || r.lo > hi) continue;
        const c = (X(Math.min(r.hi, hi)) + X(Math.max(r.lo, lo))) / 2, w = intH(r).toFixed(2).length * fs * 0.6;
        let row = 0; while (row < 2 && c - w / 2 < ends[row] + fs * 0.3) row++;
        ends[row] = c + w / 2; nrows = Math.max(nrows, row + 1); place.push({ r, c, row });
      }
    }
    const bandInt = showInt ? fs * (1.3 + 1.1 * nrows) : 0, titleH = fs * 1.7;
    const axisY = H - (tickH + bandInt + titleH);
    const top0 = (o.title ? fs * 2.2 : fs * 0.8) + (o.info ? fs * 1.3 : 0);
    const exps = S.exps.filter(e => e.hi > e.lo).sort((a, b) => b.hi - a.hi);
    const expH = exps.length ? Math.min((axisY - top0) * 0.4, fs * 22) : 0;
    const top = top0 + (expH ? expH + fs * 1.0 : 0) + (o.labels ? fs * 4.2 : fs * 0.6);
    const base = axisY - lw * 3, ymax = rangeMax(S.ppm, S.y, hi, lo) * 1.05 / S.vz, Yv = a => base - a / ymax * (base - top);
    /* trace */
    prims.push({ t: 'path', pts: tracePts(hi, lo, X, Yv, top - fs * 0.5, Math.max(2, Math.floor(xr - xl)), xl), lw, c: o.col });
    /* integral curves */
    if (showInt && o.curve && S.scale) {
      const maxH = Math.max.apply(null, S.ints.map(intH)) || 1, unit = (base - top) * 0.3 / maxH, p = S.ppm, y = S.y;
      for (const r of S.ints) {
        if (r.hi < lo || r.lo > hi) continue;
        let acc = 0; const pts = [], a = bsearch(p, r.hi), b = bsearch(p, r.lo), tot = intVal(r), hh = intH(r);
        if (!tot) continue;
        for (let i = a; i <= b && i < p.length; i++) { if (i > 0) acc += y[i] * Math.abs(p[i] - p[i - 1]); pts.push([X(p[i]), Math.max(top, base - lw * 6 - (acc / tot) * hh * unit)]); }
        prims.push({ t: 'path', pts, lw: lw * 0.8, c: o.col });
      }
    }
    /* peak labels */
    if (o.labels) {
      const cand = S.peaks.filter(q => q.ppm <= hi && q.ppm >= lo).sort((a, b) => b.h - a.h), kept = [];
      for (const q of cand) { const x = X(q.ppm); if (kept.length < 60 && kept.every(k => Math.abs(k.x - x) > fs * 1.15)) kept.push({ x, q }); }
      for (const k of kept) { const yy = Yv(k.q.h); if (yy < top - fs * 3) continue; prims.push({ t: 'text', x: k.x + fs * 0.3, y: Math.max(yy, top - fs * 0.2) - fs * 0.5, size: fs * 0.85, rot: -90, runs: [{ s: k.q.ppm.toFixed(2), k: 'n' }], c: o.fg }); }
    }
    /* axis */
    prims.push({ t: 'line', x1: xl, y1: axisY, x2: xr, y2: axisY, lw, c: o.fg });
    axisTicks(prims, X, hi, lo, axisY, fs, lw, o.fg, Math.max(4, Math.min(14, (xr - xl) / (fs * 5))));
    /* integral values */
    if (showInt) {
      const by = axisY + tickH + fs * 0.1;
      for (const q of place) {
        const r = q.r, xa = X(Math.min(r.hi, hi)), xb = X(Math.max(r.lo, lo));
        prims.push({ t: 'path', pts: [[xa, by - fs * 0.15], [xa, by], [xb, by], [xb, by - fs * 0.15]], lw: lw * 0.8, c: o.fg });
        prims.push({ t: 'text', x: q.c, y: by + fs * (1.15 + 1.1 * q.row), size: fs, anchor: 'middle', runs: [{ s: intH(r).toFixed(2), k: 'n' }], c: o.fg });
      }
    }
    prims.push({ t: 'text', x: (xl + xr) / 2, y: H - fs * 0.35, size: fs * 1.05, anchor: 'middle', runs: SC.rich('δ (ppm)'), c: o.fg });
    /* expansions: insets above the spectrum, joined to their regions by dashed lines */
    if (exps.length) {
      const n = exps.length, gap = fs * 3.6, slotW = (xr - xl) / n, eBot = top0 + expH, eAx = eBot - fs * 2.0;
      exps.forEach((e, i) => {
        const x0 = xl + i * slotW + (i ? gap / 2 : 0), x1 = xl + (i + 1) * slotW - (i < n - 1 ? gap / 2 : 0), w = x1 - x0;
        const IX = p => x0 + (e.hi - p) / (e.hi - e.lo) * w, em = rangeMax(S.ppm, S.y, e.hi, e.lo) * 1.1, IY = a => eAx - lw * 2 - a / em * (eAx - top0 - lw * 2);
        prims.push({ t: 'path', pts: tracePts(e.hi, e.lo, IX, IY, top0 - lw, Math.max(2, Math.floor(w)), x0), lw, c: o.col });
        prims.push({ t: 'line', x1: x0, y1: eAx, x2: x1, y2: eAx, lw, c: o.fg });
        axisTicks(prims, IX, e.hi, e.lo, eAx, fs * 0.85, lw, o.fg, Math.max(2, Math.min(6, w / (fs * 4.5))));
        const ya = eBot + fs * 0.1, xa = X(e.hi), xb = X(e.lo);
        prims.push({ t: 'line', x1: x0, y1: ya, x2: xa, y2: top - fs * 0.2, lw: lw * 0.7, c: o.fg, dash: [3, 3] });
        prims.push({ t: 'line', x1: x1, y1: ya, x2: xb, y2: top - fs * 0.2, lw: lw * 0.7, c: o.fg, dash: [3, 3] });
        prims.push({ t: 'path', pts: [[xa, top - fs * 0.45], [xa, top - fs * 0.2], [xb, top - fs * 0.2], [xb, top - fs * 0.45]], lw: lw * 0.7, c: o.fg });
      });
    }
    /* heading */
    let ty = fs * 1.3;
    if (o.title) { prims.push({ t: 'text', x: W / 2, y: ty, size: fs * 1.15, bold: true, anchor: 'middle', runs: SC.rich(o.title), c: o.fg }); ty += fs * 1.4; }
    if (o.info) prims.push({ t: 'text', x: xr, y: ty, size: fs, anchor: 'end', runs: SC.rich(headLine()), c: o.fg });
    return { prims, lay: { xl, xr, top, axisY, hi, lo, X, W, H, top0, expH } };
  }

  function draw() {
    if (!S.y) return;
    const W = stage.clientWidth, H = stage.clientHeight - 22; if (W < 50 || H < 50) return;
    const small = W < 520, font = css('--ui');
    mw.font = font;
    const b = build(W, H, { fs: small ? 10 : 11.5, lw: 1.1, col: css('--trace'), fg: css('--muted'), labels: $('e-peaks').checked, ints: $('e-ints').checked, curve: $('e-curve').checked, info: $('e-info').checked, title: $('etitle').value });
    lay = b.lay;
    if (svgEl) svgEl.remove();
    const t = document.createElement('div'); t.innerHTML = SC.toSVG(b.prims, W, H, { font });
    svgEl = t.firstChild; svgEl.setAttribute('preserveAspectRatio', 'none'); svgEl.style.height = H + 'px'; svgEl.style.top = '0';
    stage.insertBefore(svgEl, stage.firstChild);
  }
  const redraw = () => { cancelAnimationFrame(raf); raf = requestAnimationFrame(draw); };

  /* ---------- interaction ---------- */
  let drag = null, sel = null;
  const pxToPpm = x => lay.hi - (x - lay.xl) / (lay.xr - lay.xl) * (lay.hi - lay.lo);
  function selBox(a, b, mode) {
    if (!sel) { sel = document.createElement('div'); sel.style.cssText = 'position:absolute;pointer-events:none;top:0;border-left:1px solid;border-right:1px solid'; stage.appendChild(sel); }
    const c = mode === 'int' ? ['rgba(108,203,138,.18)', 'var(--good)'] : mode === 'mask' ? ['rgba(242,124,113,.18)', 'var(--danger)'] : mode === 'exp' ? ['rgba(230,163,61,.18)', 'var(--warn)'] : ['var(--accent-soft)', 'var(--accent)'];
    sel.style.background = c[0]; sel.style.borderColor = c[1];
    sel.style.left = Math.min(a, b) + 'px'; sel.style.width = Math.abs(b - a) + 'px'; sel.style.height = (lay.axisY + 2) + 'px'; sel.hidden = false;
  }
  function setMode(m) {
    S.mode = S.mode === m ? null : m;
    $('b-addexp').classList.toggle('on', S.mode === 'exp'); $('b-addmask').classList.toggle('on', S.mode === 'mask');
    if (S.mode) toast(S.mode === 'exp' ? 'Drag across the region to expand' : 'Drag across the region to suppress');
  }
  stage.addEventListener('pointerdown', e => {
    if (!lay || !S.y) return;
    e.preventDefault(); getSelection().removeAllRanges();
    const r = stage.getBoundingClientRect(), x = e.clientX - r.left;
    if (e.button === 1) { drag = { mode: 'pan', x, view: S.view.slice() }; stage.setPointerCapture(e.pointerId); return; }
    if (e.button !== 0) return;
    drag = { mode: e.shiftKey ? 'int' : e.altKey ? 'exp' : S.mode || 'zoom', x, x0: x, moved: false }; stage.setPointerCapture(e.pointerId);
  });
  stage.addEventListener('pointermove', e => {
    if (!lay || !S.y) return;
    const r = stage.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    if (drag) {
      if (drag.mode === 'pan') { const dppm = (x - drag.x) / (lay.xr - lay.xl) * (drag.view[0] - drag.view[1]); S.view = [drag.view[0] + dppm, drag.view[1] + dppm]; redraw(); return; }
      if (Math.abs(x - drag.x0) > 3) drag.moved = true;
      if (drag.moved) selBox(drag.x0, x, drag.mode);
      return;
    }
    const tip = $('tip');
    if (x > lay.xl && x < lay.xr) { tip.hidden = false; tip.textContent = pxToPpm(x).toFixed(3) + ' ppm' + (S.refMode ? ' · click a peak' : ''); tip.style.left = Math.min(x + 12, r.width - 150) + 'px'; tip.style.top = Math.max(4, y - 26) + 'px'; } else tip.hidden = true;
  });
  stage.addEventListener('pointerleave', () => { $('tip').hidden = true; });
  stage.addEventListener('pointerup', e => {
    if (!drag || !lay) return;
    const r = stage.getBoundingClientRect(), x = e.clientX - r.left, d = drag; drag = null;
    if (sel) sel.hidden = true;
    if (d.mode === 'pan') return;
    if (!d.moved) { if (S.refMode) clickRef(x); return; }
    const a = pxToPpm(Math.max(lay.xl, Math.min(lay.xr, d.x0))), b = pxToPpm(Math.max(lay.xl, Math.min(lay.xr, x)));
    const hi = Math.max(a, b), lo = Math.min(a, b);
    if (hi - lo < 1e-4) return;
    if (d.mode === 'int') addInt(lo, hi);
    else if (d.mode === 'exp') { if (S.exps.length >= 4) toast('Up to four expansions. Remove one first.'); else { S.exps.push({ hi, lo }); renderChips(); draw(); } if (S.mode === 'exp') setMode('exp'); }
    else if (d.mode === 'mask') { S.masks.push({ hi, lo }); if (S.mode === 'mask') setMode('mask'); finishSpectrum(false); }
    else { S.view = [hi, lo]; draw(); }
  });
  stage.addEventListener('dblclick', () => { if (S.y) { S.view = fullView(); S.vz = 1; draw(); } });
  stage.addEventListener('selectstart', e => e.preventDefault());
  stage.addEventListener('wheel', e => {
    if (!S.y) return; e.preventDefault();
    if (e.ctrlKey || e.metaKey) { const sp = S.view[0] - S.view[1], d = (e.deltaY > 0 ? 1 : -1) * sp * 0.08; S.view = [S.view[0] + d, S.view[1] + d]; }
    else S.vz = Math.max(0.05, Math.min(400, S.vz * Math.exp(-e.deltaY * 0.0015)));
    redraw();
  }, { passive: false });
  window.addEventListener('resize', redraw);
  new ResizeObserver(redraw).observe(stage);

  function clickRef(x) {
    const ppm = pxToPpm(x), tol = 8 / (lay.xr - lay.xl) * (lay.hi - lay.lo);
    let best = null; for (const p of S.rawPeaks) if (Math.abs(p.ppm - ppm) < tol && (!best || p.h > best.h)) best = p;
    if (!best) return toast('No peak there. Click closer to a peak.');
    const to = parseFloat($('refto').value) || 0;
    S.refShift += to - best.ppm; S.refMode = false; $('b-ref').classList.remove('on'); $('b-refclr').hidden = false;
    finishSpectrum(false); toast('Peak set to ' + to + ' ppm');
  }

  /* ---------- resizable panels ---------- */
  function splitter(el, axis, target, prop, min, maxFn, dflt) {
    let st = null;
    const cur = () => parseFloat(getComputedStyle(target).getPropertyValue(prop)) || dflt;
    const set = v => { v = Math.max(min, Math.min(maxFn(), v)); target.style.setProperty(prop, v + 'px'); redraw(); };
    el.addEventListener('pointerdown', e => { e.preventDefault(); st = { p: axis === 'x' ? e.clientX : e.clientY, v: cur() }; el.setPointerCapture(e.pointerId); el.classList.add('drag'); document.body.classList.add('dragging'); });
    el.addEventListener('pointermove', e => { if (!st) return; const d = (axis === 'x' ? e.clientX : e.clientY) - st.p; set(axis === 'x' ? st.v + d : st.v - d); });
    const end = () => { st = null; el.classList.remove('drag'); document.body.classList.remove('dragging'); };
    el.addEventListener('pointerup', end); el.addEventListener('pointercancel', end);
    el.addEventListener('dblclick', () => { target.style.removeProperty(prop); redraw(); });
    el.addEventListener('keydown', e => {
      const dx = { ArrowLeft: -20, ArrowRight: 20 }[e.key], dy = { ArrowUp: 20, ArrowDown: -20 }[e.key], d = axis === 'x' ? dx : dy;
      if (d === undefined) return; e.preventDefault(); set(cur() + d);
    });
  }
  splitter($('sp-v'), 'x', $('main'), '--sw', 230, () => Math.min(620, $('main').clientWidth - 360), 310);
  splitter($('sp-h'), 'y', $('stagecol'), '--lh', 80, () => $('stagecol').clientHeight - 200, 230);

  /* ---------- integrals ---------- */
  function addInt(lo, hi) {
    const r = { lo, hi }; S.ints.push(r); S.ints.sort((a, b) => b.hi - a.hi);
    if (!S.scale) S.scale = intVal(r) || 1;
    renderTable(); draw();
  }
  function autoInt() {
    if (!S.y) return;
    const p = S.ppm, y = S.y, th = Math.max(S.noise * 10, Math.max.apply(null, y) * 0.004), gap = S.nuc === '1H' ? 0.06 : 0.5, pad = S.nuc === '1H' ? 0.02 : 0.1, regs = [];
    let cur = null;
    for (let i = 0; i < y.length; i++) {
      if (y[i] > th) { if (cur && cur.lo - p[i] < gap) cur.lo = p[i]; else { cur = { hi: p[i], lo: p[i] }; regs.push(cur); } }
    }
    let out = regs.map(r => ({ hi: r.hi + pad, lo: r.lo - pad }));
    const vals = out.map(intVal), mx = Math.max.apply(null, vals) || 1;
    out = out.filter((r, i) => vals[i] > mx * 0.004);
    S.ints = out.sort((a, b) => b.hi - a.hi);
    S.scale = null;
    if (S.ints.length) { const vs = S.ints.map(intVal); S.scale = Math.min.apply(null, vs.filter(v => v > mx * 0.1).concat([mx])) || 1; }
    renderTable(); draw();
  }
  function renderTable() {
    $('inone').hidden = S.ints.length > 0;
    $('itb').innerHTML = S.ints.map((r, i) => '<tr><td class="n">' + r.hi.toFixed(2) + ' – ' + r.lo.toFixed(2) + '</td><td class="n"><input type="number" step="0.1" min="0" data-i="' + i + '" value="' + intH(r).toFixed(2) + '"> H</td><td><button class="x" data-d="' + i + '" title="Remove">✕</button></td></tr>').join('');
    report();
  }
  $('itb').addEventListener('change', e => { const i = e.target.dataset.i; if (i === undefined) return; const v = parseFloat(e.target.value); if (v > 0) { S.scale = intVal(S.ints[+i]) / v; renderTable(); draw(); } });
  $('itb').addEventListener('click', e => { const d = e.target.dataset.d; if (d === undefined) return; S.ints.splice(+d, 1); if (!S.ints.length) S.scale = null; renderTable(); draw(); });
  $('chips').addEventListener('click', e => {
    const x = e.target.dataset.e, m = e.target.dataset.m;
    if (x !== undefined) { S.exps.splice(+x, 1); renderChips(); draw(); }
    if (m !== undefined) { S.masks.splice(+m, 1); finishSpectrum(false); }
  });

  /* ---------- report (ACS style) ---------- */
  function report() {
    const el = $('acs');
    if (!S.y) { el.textContent = '—'; return; }
    const m = mhz(), sol = (S.fid ? S.fid.solvent : S.meta.solvent) || '', head = nucLabel(S.nuc) + ' NMR (' + (m ? m + ' MHz' : '') + (m && sol ? ', ' : '') + solvLabel(sol.replace(/^DMSO$/i, 'DMSO-d6')) + ') δ ';
    let items;
    if (S.nuc === '1H' && S.ints.length) {
      items = S.ints.map(r => {
        const pk = S.peaks.filter(q => q.ppm <= r.hi && q.ppm >= r.lo), h = intH(r), hs = Math.abs(h - Math.round(h)) < 0.15 ? String(Math.round(h)) : h.toFixed(1);
        const pos = pk.length <= 1 ? (pk[0] ? pk[0].ppm : (r.hi + r.lo) / 2).toFixed(2) : r.hi.toFixed(2) + '–' + r.lo.toFixed(2).replace(/^-/, '−');
        return pos + ' (' + (pk.length <= 1 ? 's' : 'm') + ', ' + hs + 'H)';
      });
    } else items = S.peaks.slice().sort((a, b) => b.ppm - a.ppm).map(q => q.ppm.toFixed(S.nuc === '1H' ? 2 : 1));
    S.report = head + items.join(', ') + '.';
    el.textContent = S.report;
  }
  $('b-copy').onclick = () => { if (!S.report) return; (navigator.clipboard ? navigator.clipboard.writeText(S.report) : Promise.reject()).then(() => toast('Copied'), () => { const r = document.createRange(); r.selectNodeContents($('acs')); const s = getSelection(); s.removeAllRanges(); s.addRange(r); toast('Press Ctrl+C to copy'); }); };

  /* ---------- export ---------- */
  function exportScene() {
    const mmW = +$('ew').value, mmH = Math.round(mmW * +$('eh').value), k = 96 / 25.4, W = Math.round(mmW * k), H = Math.round(mmH * k), pt = 96 / 72;
    const font = $('efont').value; mw.font = font;
    const b = build(W, H, { fs: +$('efs').value * pt, lw: +$('elw').value * pt, col: $('ecol').value, fg: '#000000', labels: $('e-peaks').checked, ints: $('e-ints').checked, curve: $('e-curve').checked, info: $('e-info').checked, title: $('etitle').value });
    return { prims: b.prims, W, H, mm: [mmW, mmH], font };
  }
  function save(blob, name) { const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 4000); }
  const baseName = () => (S.title || 'spectrum').replace(/\.[^.]*$/, '').replace(/[^\w.-]+/g, '_').replace(/_?fid$/i, '') || 'spectrum';
  $('b-pdf').onclick = () => {
    if (!S.y) return toast('Open data first.');
    const e = exportScene(), serif = /Times/.test(e.font);
    const u8 = SC.toPDF(e.prims, e.W, e.H, { mm: e.mm, serif, title: baseName(), bg: $('e-bg').checked ? '#ffffff' : null, mw: (s, sz, b) => mw(s, sz, b, e.font) });
    save(new Blob([u8], { type: 'application/pdf' }), baseName() + '.pdf');
  };
  $('b-svg').onclick = () => { if (!S.y) return toast('Open data first.'); const e = exportScene(); save(new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + SC.toSVG(e.prims, e.W, e.H, { font: e.font, mm: e.mm, bg: $('e-bg').checked ? '#ffffff' : null })], { type: 'image/svg+xml' }), baseName() + '.svg'); };
  $('b-png').onclick = () => {
    if (!S.y) return toast('Open data first.');
    const e = exportScene(), dpi = +$('edpi').value, cw = Math.round(e.mm[0] / 25.4 * dpi), ch = Math.round(e.mm[1] / 25.4 * dpi);
    const img = new Image();
    img.onload = () => { const c = document.createElement('canvas'); c.width = cw; c.height = ch; const g = c.getContext('2d'); if ($('e-bg').checked) { g.fillStyle = '#fff'; g.fillRect(0, 0, cw, ch); } g.drawImage(img, 0, 0, cw, ch); c.toBlob(bl => save(bl, baseName() + '_' + dpi + 'dpi.png'), 'image/png'); };
    img.onerror = () => toast('Could not render the PNG. Try PDF or SVG.');
    img.src = URL.createObjectURL(new Blob([SC.toSVG(e.prims, e.W, e.H, { font: e.font, mm: e.mm })], { type: 'image/svg+xml' }));
  };

  /* ---------- wiring ---------- */
  const reproc = () => recompute(false);
  let t1 = 0; const lazy = f => () => { cancelAnimationFrame(t1); t1 = requestAnimationFrame(f); };
  $('lb').oninput = lazy(reproc); $('zf').onchange = reproc; $('base').onchange = reproc; $('bord').onchange = reproc; $('flip').onchange = () => recompute(true);
  $('p0').oninput = $('p1').oninput = lazy(() => { S.auto = false; $('b-auto').classList.remove('on'); reproc(); });
  $('b-auto').onclick = () => { S.auto = true; $('b-auto').classList.add('on'); reproc(); };
  $('snr').oninput = lazy(() => { $('snrv').textContent = $('snr').value + '×'; if (S.y) { S.rawPeaks = R.pickPeaks(S.ppm, S.yRaw, { snr: +$('snr').value }); S.peaks = R.pickPeaks(S.ppm, S.y, { snr: +$('snr').value }); report(); draw(); } });
  $('b-ref').onclick = () => { if (!S.y) return; S.refMode = !S.refMode; $('b-ref').classList.toggle('on', S.refMode); if (S.refMode) toast('Click the peak to reference'); };
  $('b-refclr').onclick = () => { S.refShift = 0; $('b-refclr').hidden = true; finishSpectrum(false); };
  $('b-solv').onclick = setSolvent;
  $('solv').onchange = () => { if (S.yRaw && ($('m-solv').checked || $('m-imp').checked)) finishSpectrum(false); };
  for (const id of ['m-solv', 'm-imp']) $(id).onchange = () => { if (!S.yRaw) return; if ($(id).checked && !$('solv').value) toast('Choose the solvent first. Reference the spectrum before suppressing.'); finishSpectrum(false); };
  $('b-addexp').onclick = () => { if (S.y) setMode('exp'); }; $('b-addmask').onclick = () => { if (S.y) setMode('mask'); };
  $('b-autoint').onclick = autoInt; $('b-clrint').onclick = () => { S.ints = []; S.scale = null; renderTable(); draw(); };
  $('b-addint').onclick = () => { if (S.view) { const mid = (S.view[0] + S.view[1]) / 2, w = (S.view[0] - S.view[1]) * 0.1; addInt(mid - w, mid + w); toast('Added in the middle of the view. Shift+drag to place ranges directly.'); } };
  $('elw').oninput = () => { $('elwv').textContent = (+$('elw').value).toFixed(2) + ' pt'; };
  for (const id of ['e-peaks', 'e-ints', 'e-curve', 'e-info']) $(id).onchange = draw;
  $('etitle').oninput = draw;
  $('b-files').onclick = () => $('f-files').click(); $('b-folder').onclick = () => $('f-folder').click();
  $('f-files').onchange = e => { load([...e.target.files].map(f => ({ path: f.name, file: f }))); e.target.value = ''; };
  $('f-folder').onchange = e => { load([...e.target.files].map(f => ({ path: f.webkitRelativePath || f.name, file: f }))); e.target.value = ''; };
  const dz = $('drop');
  ['dragenter', 'dragover'].forEach(t => window.addEventListener(t, e => { e.preventDefault(); dz.classList.add('over'); }));
  ['dragleave', 'drop'].forEach(t => window.addEventListener(t, e => { if (t === 'drop' || e.target === document.documentElement) dz.classList.remove('over'); }));
  window.addEventListener('drop', async e => { e.preventDefault(); dz.classList.remove('over'); load(await fromDrop(e.dataTransfer)); });
  new MutationObserver(redraw).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  setControls(false); fillSolvents(); renderChips();

  /* test hook */
  window.NMRReadApp = { S, load, recompute, build, exportScene, autoInt, draw, finishSpectrum };
})();
