/* CCX NMRRead - page logic: load, process, view, integrate, export. */
(function () {
  'use strict';
  const R = window.NMRRead, $ = id => document.getElementById(id);
  const S = { fid: null, ppm0: null, y: null, ppm: null, view: null, vz: 1, ints: [], scale: null, refShift: 0, peaks: [], refMode: false, nuc: '1H', meta: {}, kind: '' };
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
      S.fid = fid; S.kind = 'fid'; S.title = title; S.refShift = 0; S.ints = []; S.scale = null; S.vz = 1;
      $('b-auto').classList.add('on'); S.auto = true; $('flip').checked = false;
      S.nuc = fid.nuc;
      fillSolvents(); preselectSolvent(fid.solvent);
      setControls(true);
      recompute(true);
    } catch (e) { toast(e.message || String(e)); console.error(e); }
  }
  function setSpectrum(d, title) {
    const sp = R.fromSpectrum(d);
    S.fid = null; S.kind = 'spectrum'; S.title = title; S.ppm0 = sp.ppm; S.y = sp.y; S.nuc = d.nuc; S.refShift = 0; S.ints = []; S.scale = null; S.vz = 1;
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
    S.ppm0 = r.ppm; S.y = r.y;
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
    S.noise = R.noiseLevel(S.y);
    S.peaks = R.pickPeaks(S.ppm, S.y, { snr: +$('snr').value });
    if (resetView || !S.view) S.view = fullView();
    $('msg').hidden = true;
    renderInfo(); renderTable(); draw();
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
    const win = S.nuc === '13C' ? 1.0 : 0.35, c = S.peaks.filter(p => Math.abs(p.ppm - t) < win);
    if (!c.length) return toast('No peak found near the solvent signal. Lower the peak threshold or use Pick peak.');
    const mx = Math.max.apply(null, c.map(p => p.h)), use = c.filter(p => p.h > mx * 0.3);
    const cen = use.reduce((a, p) => a + p.ppm * p.h, 0) / use.reduce((a, p) => a + p.h, 0);
    S.refShift += t - cen; finishSpectrum(false); toast('Referenced to ' + k + ' at ' + t + ' ppm');
    $('b-refclr').hidden = false;
  }

  /* ---------- scene ---------- */
  function niceStep(span, target) {
    const raw = span / target, pw = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / pw;
    const f = m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10;
    return { step: f * pw, minor: f === 2 ? 4 : 5 };
  }
  function viewMax(v) {
    let m = 0; const p = S.ppm, y = S.y;
    for (let i = 0; i < p.length; i++) if (p[i] <= v[0] && p[i] >= v[1] && y[i] > m) m = y[i];
    return m || 1;
  }
  function intVal(r) { return R.integrate(S.ppm, S.y, r.lo, r.hi); }
  function intH(r) { return S.scale ? intVal(r) / S.scale : 0; }

  /* Build the spectrum as an SVG string. o: {fs, lw, font, col, fg, bg, labels, ints, curve, info, title} */
  function build(W, H, o) {
    const v = S.view, hi = v[0], lo = v[1], fs = o.fs, lw = o.lw;
    const showInt = o.ints && S.ints.length, tickH = fs * 1.9;
    const xl = fs * 1.2, xr = W - fs * 1.2, span = hi - lo;
    const X = p => xl + (hi - p) / span * (xr - xl);
    /* integral labels: place each on the first row where it does not touch the previous one */
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
    const top = (o.title ? fs * 2.2 : fs * 0.8) + (o.info ? fs * 1.3 : 0) + (o.labels ? fs * 4.2 : fs * 0.6);
    const axisY = H - (tickH + bandInt + titleH);
    const base = axisY - lw * 3, ymax = viewMax(v) * 1.05 / S.vz, Yv = a => base - a / ymax * (base - top);
    const parts = [];
    if (o.bg) parts.push('<rect width="' + W + '" height="' + H + '" fill="' + o.bg + '"/>');
    parts.push('<defs><clipPath id="clp"><rect x="' + xl + '" y="' + (top - fs * 0.5) + '" width="' + (xr - xl) + '" height="' + (axisY - top + fs * 0.5) + '"/></clipPath></defs>');
    /* trace */
    const p = S.ppm, y = S.y, n = p.length;
    const i0 = Math.max(0, bsearch(p, hi) - 2), i1 = Math.min(n - 1, bsearch(p, lo) + 2), cols = Math.max(2, Math.floor(xr - xl));
    let d = '';
    if (i1 - i0 < cols * 3) { for (let i = i0; i <= i1; i++) d += (i === i0 ? 'M' : 'L') + X(p[i]).toFixed(2) + ' ' + Yv(y[i]).toFixed(2); }
    else {
      const per = (i1 - i0 + 1) / cols;
      for (let c = 0; c < cols; c++) {
        const a = Math.floor(i0 + c * per), b = Math.min(i1, Math.floor(i0 + (c + 1) * per));
        let mn = Infinity, mx = -Infinity; for (let i = a; i <= b; i++) { if (y[i] < mn) mn = y[i]; if (y[i] > mx) mx = y[i]; }
        const x = (xl + c + 0.5).toFixed(2), up = c % 2 === 0;
        d += (c === 0 ? 'M' : 'L') + x + ' ' + Yv(up ? mn : mx).toFixed(2) + 'L' + x + ' ' + Yv(up ? mx : mn).toFixed(2);
      }
    }
    parts.push('<g clip-path="url(#clp)"><path d="' + d + '" fill="none" stroke="' + o.col + '" stroke-width="' + lw + '" stroke-linejoin="round"/></g>');
    /* integral curves */
    if (showInt && o.curve && S.scale) {
      const maxH = Math.max.apply(null, S.ints.map(intH)) || 1, unit = (base - top) * 0.3 / maxH;
      for (const r of S.ints) {
        if (r.hi < lo || r.lo > hi) continue;
        let acc = 0, dd = '', first = true; const a = bsearch(p, r.hi), b = bsearch(p, r.lo), tot = intVal(r), hh = intH(r);
        if (!tot) continue;
        for (let i = a; i <= b && i < n; i++) { if (i > 0) acc += y[i] * Math.abs(p[i] - p[i - 1]); const yy = base - lw * 6 - (acc / tot) * hh * unit; dd += (first ? 'M' : 'L') + X(p[i]).toFixed(2) + ' ' + yy.toFixed(2); first = false; }
        parts.push('<path d="' + dd + '" fill="none" stroke="' + o.col + '" stroke-width="' + lw * 0.8 + '" opacity=".75" clip-path="url(#clp)"/>');
      }
    }
    /* peak labels */
    if (o.labels) {
      const cand = S.peaks.filter(q => q.ppm <= hi && q.ppm >= lo).sort((a, b) => b.h - a.h), kept = [];
      for (const q of cand) { const x = X(q.ppm); if (kept.length < 60 && kept.every(k => Math.abs(k.x - x) > fs * 1.15)) kept.push({ x, q }); }
      for (const k of kept) { const yy = Yv(k.q.h); if (yy < top - fs * 3) continue; parts.push('<text transform="translate(' + (k.x + fs * 0.35).toFixed(2) + ' ' + (Math.max(yy, top - fs * 0.2) - fs * 0.5).toFixed(2) + ') rotate(-90)" font-size="' + (fs * 0.85).toFixed(2) + '" fill="' + o.fg + '">' + k.q.ppm.toFixed(2) + '</text>'); }
    }
    /* axis */
    parts.push('<line x1="' + xl + '" x2="' + xr + '" y1="' + axisY + '" y2="' + axisY + '" stroke="' + o.fg + '" stroke-width="' + lw + '"/>');
    const nt = niceStep(span, Math.max(4, Math.min(14, (xr - xl) / (fs * 5))));
    const minor = nt.step / nt.minor, dec = Math.max(0, -Math.floor(Math.log10(nt.step) + 1e-9));
    for (let t = Math.ceil(lo / minor - 1e-9) * minor; t <= hi + 1e-9; t += minor) {
      const major = Math.abs(t / nt.step - Math.round(t / nt.step)) < 1e-6, x = X(t);
      parts.push('<line x1="' + x.toFixed(2) + '" x2="' + x.toFixed(2) + '" y1="' + axisY + '" y2="' + (axisY + (major ? fs * 0.55 : fs * 0.3)) + '" stroke="' + o.fg + '" stroke-width="' + lw + '"/>');
      if (major) parts.push('<text x="' + x.toFixed(2) + '" y="' + (axisY + fs * 1.65).toFixed(2) + '" font-size="' + fs + '" text-anchor="middle" fill="' + o.fg + '">' + (Math.abs(t) < 1e-9 ? 0 : t).toFixed(dec) + '</text>');
    }
    /* integral values */
    if (showInt) {
      const by = axisY + tickH + fs * 0.1;
      for (const q of place) {
        const r = q.r, xa = X(Math.min(r.hi, hi)), xb = X(Math.max(r.lo, lo));
        parts.push('<path d="M' + xa.toFixed(2) + ' ' + (by - fs * 0.15) + 'V' + by + 'H' + xb.toFixed(2) + 'V' + (by - fs * 0.15) + '" fill="none" stroke="' + o.fg + '" stroke-width="' + lw * 0.8 + '"/>');
        parts.push('<text x="' + q.c.toFixed(2) + '" y="' + (by + fs * (1.15 + 1.1 * q.row)).toFixed(2) + '" font-size="' + fs + '" text-anchor="middle" fill="' + o.fg + '">' + intH(r).toFixed(2) + '</text>');
      }
    }
    parts.push('<text x="' + ((xl + xr) / 2).toFixed(2) + '" y="' + (H - fs * 0.35).toFixed(2) + '" font-size="' + (fs * 1.05).toFixed(2) + '" text-anchor="middle" fill="' + o.fg + '">δ (ppm)</text>');
    /* heading */
    let ty = fs * 1.3;
    if (o.title) { parts.push('<text x="' + (W / 2) + '" y="' + ty.toFixed(2) + '" font-size="' + (fs * 1.15).toFixed(2) + '" font-weight="600" text-anchor="middle" fill="' + o.fg + '">' + esc(o.title) + '</text>'); ty += fs * 1.4; }
    if (o.info) { const t = infoLine(); if (t) parts.push('<text x="' + (xr).toFixed(2) + '" y="' + ty.toFixed(2) + '" font-size="' + fs + '" text-anchor="end" fill="' + o.fg + '">' + esc(t) + '</text>'); }
    const svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' + W + ' ' + H + '" font-family="' + o.font + '"' + (o.mm ? ' width="' + o.mm[0] + 'mm" height="' + o.mm[1] + 'mm"' : '') + '>' + parts.join('') + '</svg>';
    return { svg, lay: { xl, xr, top, axisY, hi, lo, X, Yv, W, H } };
  }
  function bsearch(p, val) { /* first index with p[i] <= val (p descending) */
    let a = 0, b = p.length - 1;
    while (a < b) { const m = (a + b) >> 1; if (p[m] <= val) b = m; else a = m + 1; }
    return a;
  }
  function infoLine() {
    const m = mhz(), sol = (S.fid ? S.fid.solvent : S.meta.solvent) || '';
    return nucLabel(S.nuc) + ' NMR' + (m || sol ? ' (' + (m ? m + ' MHz' : '') + (m && sol ? ', ' : '') + solvLabel(sol.replace(/^DMSO$/i, 'DMSO-d6')) + ')' : '');
  }

  function draw() {
    if (!S.y) return;
    const W = stage.clientWidth, H = stage.clientHeight; if (!W || !H) return;
    const small = W < 520;
    const b = build(W, H - 16, { fs: small ? 10 : 11.5, lw: 1.1, font: css('--ui'), col: css('--trace'), fg: css('--muted'), bg: '', labels: $('e-peaks').checked, ints: $('e-ints').checked, curve: $('e-curve').checked, info: $('e-info').checked, title: $('etitle').value });
    lay = b.lay;
    if (svgEl) svgEl.remove();
    const t = document.createElement('div'); t.innerHTML = b.svg; svgEl = t.firstChild; svgEl.style.fontFamily = css('--ui');
    svgEl.setAttribute('preserveAspectRatio', 'none'); stage.insertBefore(svgEl, stage.firstChild);
    svgEl.style.height = (H - 16) + 'px'; svgEl.style.top = '0';
  }
  const redraw = () => { cancelAnimationFrame(raf); raf = requestAnimationFrame(draw); };

  /* ---------- interaction ---------- */
  let drag = null, sel = null;
  const pxToPpm = x => lay.hi - (x - lay.xl) / (lay.xr - lay.xl) * (lay.hi - lay.lo);
  function selBox(a, b, cls) {
    if (!sel) { sel = document.createElement('div'); sel.style.cssText = 'position:absolute;pointer-events:none;top:0;border-left:1px solid var(--accent);border-right:1px solid var(--accent)'; stage.appendChild(sel); }
    sel.style.background = cls === 'int' ? 'rgba(108,203,138,.18)' : 'var(--accent-soft)'; sel.style.borderColor = cls === 'int' ? 'var(--good)' : 'var(--accent)';
    sel.style.left = Math.min(a, b) + 'px'; sel.style.width = Math.abs(b - a) + 'px'; sel.style.height = (lay.axisY + 2) + 'px'; sel.hidden = false;
  }
  stage.addEventListener('pointerdown', e => {
    if (!lay || !S.y) return;
    const r = stage.getBoundingClientRect(), x = e.clientX - r.left;
    if (e.button === 1) { drag = { mode: 'pan', x, view: S.view.slice() }; stage.setPointerCapture(e.pointerId); e.preventDefault(); return; }
    if (e.button !== 0) return;
    drag = { mode: e.shiftKey ? 'int' : 'zoom', x, x0: x, moved: false }; stage.setPointerCapture(e.pointerId);
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
    if (d.mode === 'int') addInt(lo, hi); else { S.view = [hi, lo]; draw(); }
  });
  stage.addEventListener('dblclick', () => { if (S.y) { S.view = fullView(); S.vz = 1; draw(); } });
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
    let best = null; for (const p of S.peaks) if (Math.abs(p.ppm - ppm) < tol && (!best || p.h > best.h)) best = p;
    if (!best) return toast('No peak there. Click closer to a peak.');
    const to = parseFloat($('refto').value) || 0;
    S.refShift += to - best.ppm; S.refMode = false; $('b-ref').classList.remove('on'); $('b-refclr').hidden = false;
    finishSpectrum(false); toast('Peak set to ' + to + ' ppm');
  }

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
    return build(W, H, { fs: +$('efs').value * pt, lw: +$('elw').value * pt, font: $('efont').value, col: $('ecol').value, fg: '#000000', bg: '', labels: $('e-peaks').checked, ints: $('e-ints').checked, curve: $('e-curve').checked, info: $('e-info').checked, title: $('etitle').value, mm: [mmW, mmH] });
  }
  function save(blob, name) { const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 4000); }
  const baseName = () => (S.title || 'spectrum').replace(/\.[^.]*$/, '').replace(/[^\w.-]+/g, '_').replace(/_?fid$/i, '') || 'spectrum';
  $('b-svg').onclick = () => { if (!S.y) return toast('Open data first.'); const b = exportScene(), bg = $('e-bg').checked; save(new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + b.svg.replace('>', '>' + (bg ? '<rect width="100%" height="100%" fill="#fff"/>' : '')) ], { type: 'image/svg+xml' }), baseName() + '.svg'); };
  $('b-png').onclick = () => {
    if (!S.y) return toast('Open data first.');
    const b = exportScene(), dpi = +$('edpi').value, mmW = +$('ew').value, mmH = Math.round(mmW * +$('eh').value), cw = Math.round(mmW / 25.4 * dpi), ch = Math.round(mmH / 25.4 * dpi);
    const img = new Image();
    img.onload = () => { const c = document.createElement('canvas'); c.width = cw; c.height = ch; const g = c.getContext('2d'); if ($('e-bg').checked) { g.fillStyle = '#fff'; g.fillRect(0, 0, cw, ch); } g.drawImage(img, 0, 0, cw, ch); c.toBlob(bl => save(bl, baseName() + '_' + dpi + 'dpi.png'), 'image/png'); };
    img.onerror = () => toast('Could not render the PNG. Try SVG.');
    img.src = URL.createObjectURL(new Blob([b.svg], { type: 'image/svg+xml' }));
  };

  /* ---------- wiring ---------- */
  const reproc = () => recompute(false);
  let t1 = 0; const lazy = f => () => { cancelAnimationFrame(t1); t1 = requestAnimationFrame(f); };
  $('lb').oninput = lazy(reproc); $('zf').onchange = reproc; $('base').onchange = reproc; $('bord').onchange = reproc; $('flip').onchange = () => recompute(true);
  $('p0').oninput = $('p1').oninput = lazy(() => { S.auto = false; $('b-auto').classList.remove('on'); reproc(); });
  $('b-auto').onclick = () => { S.auto = true; $('b-auto').classList.add('on'); reproc(); };
  $('snr').oninput = lazy(() => { $('snrv').textContent = $('snr').value + '×'; if (S.y) { S.peaks = R.pickPeaks(S.ppm, S.y, { snr: +$('snr').value }); report(); draw(); } });
  $('b-ref').onclick = () => { if (!S.y) return; S.refMode = !S.refMode; $('b-ref').classList.toggle('on', S.refMode); if (S.refMode) toast('Click the peak to reference'); };
  $('b-refclr').onclick = () => { S.refShift = 0; $('b-refclr').hidden = true; finishSpectrum(false); };
  $('b-solv').onclick = setSolvent;
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
  $('b-theme').addEventListener('click', () => setTimeout(redraw, 50));
  new MutationObserver(redraw).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  setControls(false); fillSolvents();

  /* test hook */
  window.NMRReadApp = { S, load, recompute, build, exportScene, autoInt, draw };
})();
