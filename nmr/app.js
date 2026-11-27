/* CCX NMR page logic: embeds CCX Draw, predicts with nmr.js, draws the spectrum, table and structure. */
(function () {
'use strict';
const $ = s => document.querySelector(s);
const S = { nuc: '1H', mhz: 400, lw: 1, solv: 'CDCl3', ex: true, nums: false, res: null, err: '', lines: [], y: null, grid: null, view: null, sel: null, last: null, drag: null, hover: null, vz: 1, exps: [], armed: false, pan: null, insets: [], opt: { icurve: true, brace: true } };
const dark = () => { const t = document.documentElement.dataset.theme; return t ? t === 'dark' : matchMedia('(prefers-color-scheme:dark)').matches; };
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const toast = m => { const t = $('#toast'); t.textContent = m; t.hidden = false; clearTimeout(toast.t); toast.t = setTimeout(() => t.hidden = true, 2200); };

/* ---------------- Draw bridge (same contract as CCX Dock: global molfile()) ---------------- */
const BL = { frame: $('#bl') };
const blWin = () => { try { return BL.frame.contentWindow; } catch (e) { return null; } };
function molfile() { const w = blWin(); try { if (w && typeof w.molfile === 'function') return w.molfile(); } catch (e) {} return null; }
/* Draw inside this page: keep Draw's tall, vertical layout at its natural size and let the column scroll. */
function prepareDraw() {
  const w = blWin(); try {
    const d = w.document, st = d.createElement('style');
    st.textContent = '@media (max-width:880px){html,body{height:auto!important}#app{height:auto!important;min-height:0!important;grid-template-rows:auto auto 560px auto auto!important}}';
    d.head.appendChild(st);
    /* wheel over Draw's toolbars and panels scrolls this column; over the drawing area it still zooms the drawing */
    d.addEventListener('wheel', ev => {
      if (ev.ctrlKey || ev.defaultPrevented) return; const t = ev.target; if (t.closest && t.closest('#stage,textarea,select,.box')) return;
      const L = $('#left'); if (L.scrollHeight <= L.clientHeight) return; L.scrollTop += ev.deltaMode === 1 ? ev.deltaY * 32 : ev.deltaY; ev.preventDefault();
    }, { passive: false });
  } catch (e) {}
}
function fitDraw() {
  const L = $('#left'), w = blWin(); if (!w) return; let h = L.clientHeight;
  try {
    if (L.clientWidth > 880) h = Math.max(L.clientHeight, 640);
    else { const b = w.document.body, r = b.getBoundingClientRect(); h = Math.max(L.clientHeight, Math.ceil(r.height + parseFloat(w.getComputedStyle(b).marginTop || 0) + parseFloat(w.getComputedStyle(b).marginBottom || 0))); }
  } catch (e) {}
  if (Math.abs(BL.frame.offsetHeight - h) > 2) BL.frame.style.height = h + 'px';
}
function syncTheme() { const w = blWin(); try { if (w && w.document) w.document.documentElement.dataset.theme = dark() ? 'dark' : 'light'; } catch (e) {} }
async function loadDraw() {
  BL.frame.addEventListener('load', () => { syncTheme(); prepareDraw(); fitDraw(); poll(); });
  try { const r = await fetch('../draw/index.html'); if (!r.ok) throw new Error('HTTP ' + r.status); BL.frame.srcdoc = await r.text(); }
  catch (e) { $('#msg').textContent = 'The CCX Draw editor could not be loaded (' + e.message + ').'; }
}
function poll() { if (S.hold) return; const t = molfile(); if (t == null || t === S.last) return; S.last = t; update(); }
setInterval(() => { poll(); fitDraw(); }, 600);
window.addEventListener('cp-theme', () => { syncTheme(); draw(); drawStructure(); });

/* ---------------- prediction ---------------- */
function update() {
  S.sel = null; S.res = null; S.err = ''; S.exps = []; S.vz = 1; S.armed = false; if (typeof syncExpUI === 'function') syncExpUI();
  try { S.res = NMR.predict(S.last, { solvent: S.solv }); }
  catch (e) { S.err = /no atoms/i.test(e.message) ? '' : e.message; }
  buildSpectrum(true); renderTable(); drawStructure(); draw();
}
const full = () => S.nuc === '1H' ? [-0.5, 13.5] : [-10, 230];
function buildSpectrum(resetView) {
  const r = S.res; S.lines = []; S.y = null;
  if (!r) { $('#msg').hidden = false; $('#msg').textContent = S.err || 'Draw a structure to begin.'; return; }
  const [lo, hi] = full(); S.grid = { lo, hi, N: S.nuc === '1H' ? 24000 : 60000 };
  if (S.nuc === '1H') S.lines = NMR.hLines(r.h.signals.filter(s => !(S.solv === 'D2O' || S.solv === 'CD3OD') || !s.exch), S.mhz, { exchangeable: S.ex });
  else S.lines = NMR.cLines(r.c.signals, S.nuc === 'DEPT' ? 'dept' : 'full');
  const lw = S.nuc === '1H' ? S.lw : Math.max(S.lw, 1);
  const mhz = S.nuc === '1H' ? S.mhz : S.mhz * 0.2515;
  /* broad exchangeable lines: widen by mixing in a second pass */
  const normal = S.lines.filter(l => !l.broad), broad = S.lines.filter(l => l.broad);
  const y = NMR.synth(normal, { xmin: lo, xmax: hi, N: S.grid.N, lw, mhz });
  if (broad.length) { const yb = NMR.synth(broad, { xmin: lo, xmax: hi, N: S.grid.N, lw: lw * 8, mhz }); for (let i = 0; i < y.length; i++) y[i] += yb[i]; }
  S.y = y; S.ylw = lw; S.ymhz = mhz;
  { const cum = new Float64Array(y.length + 1); for (let i = y.length - 1; i >= 0; i--) cum[i] = cum[i + 1] + y[i]; S.cum = cum; }
  let mx = 0, mn = 0; for (let i = 0; i < y.length; i++) { if (y[i] > mx) mx = y[i]; if (y[i] < mn) mn = y[i]; }
  S.ymax = mx || 1; S.ymin = mn;
  if (resetView || !S.view) S.view = autoView();
  $('#msg').hidden = true;
}
function autoView() {
  const r = S.res; if (!r) return full();
  if (S.nuc === '1H') {
    const sh = r.h.signals.filter(s => S.ex || !s.exch).map(s => s.shift); const mx = Math.max(...sh, 8), mn = Math.min(...sh, 0);
    return [Math.max(-0.5, mn - 0.8), Math.min(13.5, mx + 0.8)];
  }
  const sh = r.c.signals.map(s => s.shift); if (!sh.length) return full();
  return [Math.max(-10, Math.min(...sh) - 12), Math.min(230, Math.max(...sh) + 12)];
}

/* ---------------- spectrum canvas ---------------- */
const cv = $('#spec'), ctx = cv.getContext('2d');
const PAD = { l: 14, r: 14, t: 26, b: 84 };
function size() { const b = $('#specbox'), d = Math.min(2, window.devicePixelRatio || 1); const w = b.clientWidth, h = b.clientHeight; if (cv.width !== w * d || cv.height !== h * d) { cv.width = w * d; cv.height = h * d; } ctx.setTransform(d, 0, 0, d, 0, 0); return [w, h]; }
const xOf = (ppm, w) => PAD.l + (S.view[1] - ppm) / (S.view[1] - S.view[0]) * (w - PAD.l - PAD.r);
const ppmOf = (px, w) => S.view[1] - (px - PAD.l) / (w - PAD.l - PAD.r) * (S.view[1] - S.view[0]);
function niceStep(range, target) { const raw = range / target, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p; return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p; }
function signalsNow() { return !S.res ? [] : S.nuc === '1H' ? S.res.h.signals : S.res.c.signals; }

const shownSigs = () => signalsNow().filter(s => !(S.nuc === '1H' && s.exch && (!S.ex || S.solv === 'D2O' || S.solv === 'CD3OD')) && !(S.nuc === 'DEPT' && s.type === 'C'));
const spanOf = s => S.nuc === '1H' && s.lines && s.lines.length ? (s.lines[s.lines.length - 1].d - s.lines[0].d) / 2 / S.mhz : 0;
/* tallest and lowest point of the simulated trace between two shifts */
function yRange(a, b) {
  const dx = (S.grid.hi - S.grid.lo) / (S.grid.N - 1), i0 = Math.max(0, Math.floor((a - S.grid.lo) / dx)), i1 = Math.min(S.grid.N - 1, Math.ceil((b - S.grid.lo) / dx));
  let mx = 0, mn = 0; for (let i = i0; i <= i1; i++) { const v = S.y[i]; if (v > mx) mx = v; if (v < mn) mn = v; }
  return { mx: mx || 1, mn, i0, i1, dx };
}
/* the spectrum trace between shifts a<b, mapped to pixels by fx */
function traceLine(fx, a, b, zero, scale, color, lw) {
  const { i0, i1, dx } = yRange(a, b), span = Math.abs(fx(a) - fx(b)), step = Math.max(1, Math.floor((i1 - i0) / (span * 2 + 1)));
  ctx.strokeStyle = color; ctx.lineWidth = lw; ctx.beginPath();
  for (let i = i0; i <= i1; i += step) { let v = S.y[i]; if (step > 1) { for (let k = 1; k < step && i + k <= i1; k++) { const q = S.y[i + k]; if (Math.abs(q) > Math.abs(v)) v = q; } } const x = fx(S.grid.lo + i * dx), y = zero - v * scale; if (i === i0) ctx.moveTo(x, y); else ctx.lineTo(x, y); }
  ctx.stroke();
}
function metrics(a, b, pTop, base, gain) {
  const { mx, mn } = yRange(a, b), ph = base - pTop, hasNeg = mn < -mx * 0.02, zero = hasNeg ? base - ph * 0.45 : base, up = hasNeg ? ph * 0.40 : ph * 0.80;
  return { zero, scale: (hasNeg ? up / Math.max(mx, -mn) : up / mx) * gain };
}
function brace(x0, x1, y) {
  const d = 5, t = 5, r = Math.max(1, Math.min(5, (x1 - x0) / 4)), m = (x0 + x1) / 2;
  ctx.beginPath(); ctx.moveTo(x0, y); ctx.quadraticCurveTo(x0, y + d, x0 + r, y + d); ctx.lineTo(m - r, y + d); ctx.quadraticCurveTo(m, y + d, m, y + d + t);
  ctx.quadraticCurveTo(m, y + d, m + r, y + d); ctx.lineTo(x1 - r, y + d); ctx.quadraticCurveTo(x1, y + d, x1, y); ctx.stroke();
}
function drawInset(e, k, ix, iy, iw, ih, c) {
  const ixl = ix + 10, ixr = ix + iw - 10, topY = iy + 22, axisY = iy + ih - 22, fx = p => ixl + (e.b - p) / (e.b - e.a) * (ixr - ixl);
  ctx.fillStyle = c.panel; ctx.strokeStyle = c.line; ctx.lineWidth = 1; ctx.beginPath(); ctx.roundRect(ix + 0.5, iy + 0.5, iw - 1, ih - 1, 8); ctx.fill(); ctx.stroke();
  const m = metrics(e.a, e.b, topY, axisY, 1);
  ctx.save(); ctx.beginPath(); ctx.rect(ixl, topY - 4, ixr - ixl, axisY - topY + 8); ctx.clip(); traceLine(fx, e.a, e.b, m.zero, m.scale, c.trace, 1.3); ctx.restore();
  ctx.strokeStyle = c.line; ctx.beginPath(); ctx.moveTo(ixl, axisY + 0.5); ctx.lineTo(ixr, axisY + 0.5); ctx.stroke();
  const st = niceStep(e.b - e.a, 4), dec = Math.max(0, Math.min(4, Math.ceil(-Math.log10(st) - 1e-9)));
  ctx.fillStyle = c.muted; ctx.font = '10.5px ' + css('--mono'); ctx.textAlign = 'center';
  for (let t = Math.ceil(e.a / st) * st; t <= e.b + 1e-9; t += st) { const x = Math.round(fx(t)) + 0.5; ctx.beginPath(); ctx.moveTo(x, axisY); ctx.lineTo(x, axisY + 4); ctx.stroke(); ctx.fillText(t.toFixed(dec), x, axisY + 16); }
  ctx.textAlign = 'left'; ctx.fillStyle = c.accent; ctx.font = '600 11px ' + css('--mono'); ctx.fillText('Expansion ' + (k + 1) + '  ' + e.b.toFixed(2) + ' – ' + e.a.toFixed(2) + ' ppm', ix + 10, iy + 15);
  ctx.textAlign = 'center'; ctx.fillStyle = c.fg; ctx.font = '11px ' + css('--mono');
  const { dx } = yRange(e.a, e.b);
  shownSigs().forEach(s => { if (s.shift < e.a || s.shift > e.b) return; const ci = Math.round((s.shift - S.grid.lo) / dx); let pk = 0; for (let q = -Math.round(0.1 / dx); q <= Math.round(0.1 / dx); q++) { const v = S.y[ci + q]; if (v !== undefined && Math.abs(v) > Math.abs(pk)) pk = v; } ctx.fillText(s.shift.toFixed(S.nuc === '1H' ? 3 : 1), fx(s.shift), Math.max(topY + 2, m.zero - pk * m.scale - 5)); });
  ctx.fillStyle = c.muted; ctx.font = '15px ' + css('--ui'); ctx.textAlign = 'center'; ctx.fillText('×', ix + iw - 14, iy + 16);
  S.insets.push({ k, x: ix, y: iy, w: iw, h: ih, cx: ix + iw - 26, cy: iy, cw: 26, ch: 24 });
}

function draw() {
  const [w, h] = size(); ctx.clearRect(0, 0, w, h); S.insets = [];
  const fg = css('--fg'), muted = css('--muted'), line = css('--line'), accent = css('--accent'), trace = css('--trace'), oneH = S.nuc === '1H';
  if (!S.y || !S.view) return;
  const nE = S.exps.length, base = h - PAD.b, EH = nE ? Math.max(90, Math.min(210, Math.min(h * 0.34, (base - PAD.t) * 0.5))) : 0, top = PAD.t + (nE ? EH + 14 : 0), ph = base - top;
  const [v0, v1] = S.view, laneOn = oneH && S.opt.icurve, laneH = laneOn ? Math.min(90, ph * 0.2) : 0, pTop = top + (laneOn ? laneH + 8 : 0);
  const mt = metrics(v0, v1, pTop, base, S.vz), zero = mt.zero, scale = mt.scale, fx = p => xOf(p, w);
  /* selection band and expansion regions */
  const sigs = signalsNow(), sel = sigs.find(s => s.id === S.sel);
  if (sel) {
    let a = sel.shift, b = sel.shift, pad = oneH ? 0.04 : 1.5;
    if (oneH) { sel.lines.forEach(l => { const p = sel.shift + l.d / S.mhz; a = Math.min(a, p); b = Math.max(b, p); }); pad = 0.05 + (sel.exch ? 0.3 : 0); }
    ctx.fillStyle = css('--sel'); const xa = xOf(b + pad, w), xb = xOf(a - pad, w); ctx.fillRect(xa, top - 4, xb - xa, base - top + 4);
  }
  /* axis */
  ctx.strokeStyle = line; ctx.fillStyle = muted; ctx.lineWidth = 1; ctx.font = '11px ' + css('--mono'); ctx.textAlign = 'center';
  ctx.beginPath(); ctx.moveTo(PAD.l, base + 0.5); ctx.lineTo(w - PAD.r, base + 0.5); ctx.stroke();
  const st = niceStep(v1 - v0, 9), minor = st / 5;
  for (let t = Math.ceil(v0 / minor) * minor; t <= v1 + 1e-9; t += minor) { const x = Math.round(xOf(t, w)) + 0.5, big = Math.abs(t / st - Math.round(t / st)) < 1e-6; ctx.beginPath(); ctx.moveTo(x, base); ctx.lineTo(x, base + (big ? 6 : 3)); ctx.stroke(); if (big) ctx.fillText((+t.toFixed(3)).toString(), x, base + 19); }
  ctx.textAlign = 'right'; ctx.fillText('δ (ppm)', w - PAD.r, 14);
  if (Math.abs(S.vz - 1) > 0.01) ctx.fillText('height ×' + (S.vz >= 10 ? S.vz.toFixed(0) : S.vz.toFixed(1)), w - PAD.r - 64, 14);
  ctx.textAlign = 'left'; ctx.fillText((oneH ? '¹H ' + S.mhz + ' MHz' : S.nuc === '13C' ? '¹³C ' + (S.mhz * 0.2515).toFixed(0) + ' MHz' : 'DEPT-135 ' + (S.mhz * 0.2515).toFixed(0) + ' MHz') + ' · ' + S.solv + ' · simulated', PAD.l, 14);
  /* trace */
  ctx.save(); ctx.beginPath(); ctx.rect(PAD.l, top - 2, w - PAD.l - PAD.r, base - top + 3); ctx.clip();
  traceLine(fx, v0, v1, zero, scale, trace, 1.3);
  if (sel) { /* redraw the selected signal in the accent colour */
    const half = oneH ? 0.06 + (sel.exch ? 0.35 : 0) + spanOf(sel) : 1.2;
    traceLine(fx, Math.max(v0, sel.shift - half), Math.min(v1, sel.shift + half), zero, scale, accent, 1.8);
  }
  ctx.restore();
  /* peak labels */
  ctx.font = '11px ' + css('--mono'); ctx.textAlign = 'center';
  const placed = [], shown = shownSigs(), dxg = (S.grid.hi - S.grid.lo) / (S.grid.N - 1);
  for (const s of shown.slice().sort((p, q) => q.shift - p.shift)) {
    if (s.shift < v0 || s.shift > v1) continue;
    const x = xOf(s.shift, w); let peak = 0; const ci = Math.round((s.shift - S.grid.lo) / dxg);
    for (let k = -Math.round(0.15 / dxg); k <= Math.round(0.15 / dxg); k++) { const v = S.y[ci + k]; if (v !== undefined && Math.abs(v) > Math.abs(peak)) peak = v; }
    const ty = peak >= 0 ? zero - peak * scale - 5 : zero - peak * scale + 13;
    if (placed.some(p => Math.abs(p - x) < 28)) continue; placed.push(x);
    ctx.fillStyle = s.id === S.sel ? accent : fg; ctx.fillText(s.shift.toFixed(oneH ? 2 : 1), x, Math.max(pTop - 2, ty));
  }
  if (oneH) {
    /* integral curve: a running integral across the top of the spectrum, one step per signal */
    if (laneOn && S.cum) {
      const inView = shown.filter(s => s.shift >= v0 && s.shift <= v1), visH = inView.reduce((t, s) => t + s.n, 0), allH = shown.reduce((t, s) => t + s.n, 0);
      if (visH > 0 && allH > 0) {
        const hY = S.cum[0] / allH, laneBase = top + laneH, pxPerH = laneH * 0.92 / visH, iEdge = Math.min(S.grid.N - 1, Math.max(0, Math.round((v1 - S.grid.lo) / dxg)));
        ctx.save(); ctx.beginPath(); ctx.rect(PAD.l, top - 4, w - PAD.l - PAD.r, laneH + 6); ctx.clip();
        ctx.strokeStyle = accent; ctx.lineWidth = 1.5; ctx.beginPath();
        for (let px = PAD.l; px <= w - PAD.r; px++) { const i = Math.min(S.grid.N - 1, Math.max(0, Math.round((ppmOf(px, w) - S.grid.lo) / dxg))), f = (S.cum[i] - S.cum[iEdge]) / hY, y = Math.max(top - 2, laneBase - f * pxPerH); if (px === PAD.l) ctx.moveTo(px, y); else ctx.lineTo(px, y); }
        ctx.stroke(); ctx.restore();
      }
    }
    /* braces and integrals under the axis */
    ctx.fillStyle = accent; ctx.strokeStyle = accent; ctx.lineWidth = 1.3; ctx.textAlign = 'center'; ctx.lineCap = 'round';
    const rowEnd = [-1e9, -1e9];
    for (const s of shown.slice().sort((p, q) => q.shift - p.shift)) {
      if (s.shift < v0 || s.shift > v1) continue;
      const half = spanOf(s) + 0.025; let x0 = Math.max(PAD.l, xOf(s.shift + half, w)), x1 = Math.min(w - PAD.r, xOf(s.shift - half, w));
      if (x1 - x0 < 16) { const c0 = (x0 + x1) / 2; x0 = c0 - 8; x1 = c0 + 8; }
      const t = s.n + 'H', tw = ctx.measureText(t).width / 2, c = (x0 + x1) / 2, left = Math.min(x0, c - tw), right = Math.max(x1, c + tw), row = left > rowEnd[0] + 6 ? 0 : 1;
      rowEnd[row] = right; const by = base + 27 + row * 22;
      if (S.opt.brace) brace(x0, x1, by);
      ctx.fillText(t, c, by + (S.opt.brace ? 22 : 11));
    }
    ctx.lineCap = 'butt';
  }
  /* expansions: regions on the spectrum, insets above, dashed connectors */
  if (nE) {
    const gap = 8, avail = w - PAD.l - PAD.r, iw = Math.min(380, (avail - gap * (nE - 1)) / nE), x00 = PAD.l + (avail - (iw * nE + gap * (nE - 1))) / 2, iy = 22;
    const col = { panel: css('--panel'), line, trace, muted, accent, fg };
    S.exps.forEach((e, k) => {
      const xa = xOf(Math.min(v1, e.b), w), xb = xOf(Math.max(v0, e.a), w), ix = x00 + k * (iw + gap);
      if (e.b >= v0 && e.a <= v1) {
        ctx.fillStyle = css('--accent-soft'); ctx.fillRect(xa, top - 4, Math.max(1, xb - xa), base - top + 4);
        ctx.strokeStyle = accent; ctx.globalAlpha = 0.55; ctx.setLineDash([4, 3]); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(ix, iy + EH); ctx.lineTo(xa, top - 4); ctx.moveTo(ix + iw, iy + EH); ctx.lineTo(xb, top - 4); ctx.stroke(); ctx.setLineDash([]); ctx.globalAlpha = 1;
        ctx.fillStyle = accent; ctx.font = '600 10.5px ' + css('--mono'); ctx.textAlign = 'center'; ctx.fillText(String(k + 1), (xa + xb) / 2, top - 8 + 0);
      }
      drawInset(e, k, ix, iy, iw, EH, col);
    });
  }
  S.geom = { w, h, zero, scale, base, top, pTop };
}
function sigAt(px) { /* nearest signal to a pixel x */
  const w = cv.clientWidth, p = ppmOf(px, w), tol = (S.view[1] - S.view[0]) * 14 / w; let best = null, bd = 1e9;
  for (const s of signalsNow()) {
    if (S.nuc === 'DEPT' && s.type === 'C') continue; if (S.nuc === '1H' && s.exch && (!S.ex || S.solv === 'D2O' || S.solv === 'CD3OD')) continue;
    const half = S.nuc === '1H' ? (s.lines[s.lines.length - 1].d - s.lines[0].d) / 2 / S.mhz : 0, d = Math.max(0, Math.abs(p - s.shift) - half);
    if (d < bd) { bd = d; best = s; }
  }
  return best && bd <= tol ? best : null;
}
const inInset = (x, y) => S.insets.find(i => x >= i.x && x <= i.x + i.w && y >= i.y && y <= i.y + i.h);
function addExp(a, b) {
  let lo = Math.min(a, b), hi = Math.max(a, b); const minW = S.nuc === '1H' ? 0.02 : 0.5;
  if (hi - lo < minW) { const c = (lo + hi) / 2; lo = c - minW / 2; hi = c + minW / 2; }
  const [fl, fh] = full(); lo = Math.max(fl, lo); hi = Math.min(fh, hi);
  if (S.exps.length >= 4) { toast('Four expansions at most. Remove one with its × first.'); return; }
  S.exps.push({ a: lo, b: hi }); S.exps.sort((p, q) => q.b - p.b); syncExpUI();
}
function syncExpUI() { $('#b-exp').classList.toggle('on', S.armed); $('#b-expclr').hidden = !S.exps.length; cv.style.cursor = S.armed ? 'cell' : 'crosshair'; }
cv.addEventListener('mousedown', e => {
  if (!S.y) return;
  if (e.button === 1) { e.preventDefault(); S.pan = { x0: e.offsetX, view: S.view.slice() }; cv.style.cursor = 'grabbing'; $('#tip').hidden = true; return; }
  if (e.button !== 0) return;
  const ins = inInset(e.offsetX, e.offsetY);
  if (ins) { if (e.offsetX >= ins.cx && e.offsetY <= ins.cy + ins.ch) { S.exps.splice(ins.k, 1); syncExpUI(); draw(); } return; }
  S.drag = { x0: e.offsetX, x1: e.offsetX, moved: false, mode: S.armed || e.altKey ? 'exp' : 'zoom' };
});
window.addEventListener('mouseup', e => {
  if (S.pan && e.button === 1) { S.pan = null; syncExpUI(); return; }
  const d = S.drag; S.drag = null; if (!d) return;
  if (d.moved && Math.abs(d.x1 - d.x0) > 6) {
    const w = cv.clientWidth, a = ppmOf(d.x0, w), b = ppmOf(d.x1, w);
    if (d.mode === 'exp') { addExp(a, b); S.armed = false; syncExpUI(); } else S.view = [Math.min(a, b), Math.max(a, b)];
    draw();
  } else if (!d.moved) { const s = sigAt(d.x0); select(s ? s.id : null); }
  $('#tip').hidden = true;
});
cv.addEventListener('auxclick', e => { if (e.button === 1) e.preventDefault(); });
cv.addEventListener('mousemove', e => {
  if (S.pan) {
    const w = cv.clientWidth, [a0, b0] = S.pan.view, range = b0 - a0, d = (e.offsetX - S.pan.x0) * range / (w - PAD.l - PAD.r), [lo, hi] = full();
    let na = a0 + d, nb = b0 + d; if (na < lo) { nb += lo - na; na = lo; } if (nb > hi) { na -= nb - hi; nb = hi; } S.view = [Math.max(lo, na), Math.min(hi, nb)]; draw(); return;
  }
  if (S.drag) {
    S.drag.x1 = e.offsetX; if (Math.abs(S.drag.x1 - S.drag.x0) > 3) S.drag.moved = true; draw();
    if (S.drag.moved) { const g = S.geom, x = Math.min(S.drag.x0, S.drag.x1), ww = Math.abs(S.drag.x1 - S.drag.x0); ctx.fillStyle = css('--accent-soft'); ctx.fillRect(x, g.top, ww, g.base - g.top); if (S.drag.mode === 'exp') { ctx.strokeStyle = css('--accent'); ctx.setLineDash([5, 3]); ctx.strokeRect(x + 0.5, g.top + 0.5, ww, g.base - g.top); ctx.setLineDash([]); } }
    return;
  }
  if (!S.y) return; const tip = $('#tip');
  if (inInset(e.offsetX, e.offsetY)) { tip.hidden = true; cv.style.cursor = e.offsetX >= S.insets[0].cx ? 'pointer' : 'default'; return; }
  const s = sigAt(e.offsetX), p = ppmOf(e.offsetX, cv.clientWidth);
  tip.hidden = false; tip.textContent = p.toFixed(S.nuc === '1H' ? 3 : 1) + ' ppm' + (s ? '  ·  ' + (S.nuc === '1H' ? s.n + 'H ' + s.mult : s.n + ' C ' + s.type) : '');
  tip.style.left = Math.min(cv.clientWidth - tip.offsetWidth - 6, e.offsetX + 12) + 'px'; tip.style.top = Math.max(4, e.offsetY - 28) + 'px'; cv.style.cursor = S.armed || e.altKey ? 'cell' : s ? 'pointer' : 'crosshair';
});
cv.addEventListener('mouseleave', () => { $('#tip').hidden = true; });
cv.addEventListener('dblclick', e => { if (inInset(e.offsetX, e.offsetY)) return; S.view = autoView(); S.vz = 1; draw(); });
/* the wheel changes peak height; a middle-button drag pans along the field axis; a plain drag zooms */
cv.addEventListener('wheel', e => {
  if (!S.y) return; e.preventDefault();
  S.vz = Math.max(0.1, Math.min(400, S.vz * (e.deltaY < 0 ? 1.2 : 1 / 1.2))); draw();
}, { passive: false });
new ResizeObserver(() => { draw(); }).observe($('#specbox'));

/* ---------------- peak table ---------------- */
function envName(s) {
  if (S.nuc !== '1H') return s.type;
  if (s.exch) return s.hs[0].env;
  return { sp3: s.label, arom: 'aromatic H', vinyl: 'alkene H', ald: 'aldehyde H', alkyne: 'alkyne H', imine: 'HC=N' }[s.kind] || s.label;
}
function Jtext(s) { return s.Jlist && s.Jlist.length ? s.Jlist.map(g => g.J.toFixed(1)).join(', ') : ''; }
function renderTable() {
  const th = $('#th'), tb = $('#tb'), note = $('#note'); tb.innerHTML = '';
  if (!S.res) { th.innerHTML = ''; note.textContent = ''; return; }
  const r = S.res, sh = (S.solv === 'D2O' || S.solv === 'CD3OD');
  if (S.nuc === '1H') {
    th.innerHTML = '<tr><th class="n">δ (ppm)</th><th class="n">H</th><th>Mult.</th><th class="n">J (Hz)</th><th>Type</th><th>Atoms</th></tr>';
    for (const s of r.h.signals) {
      if (s.exch && (!S.ex || sh)) continue;
      const tr = document.createElement('tr'); tr.className = 'row'; tr.dataset.id = s.id;
      tr.innerHTML = '<td class="n">' + s.shift.toFixed(2) + '</td><td class="n">' + s.n + '</td><td>' + s.mult + '</td><td class="n">' + Jtext(s) + '</td><td>' + envName(s) + '</td><td class="n">' + s.atoms.map(a => a + 1).join(', ') + '</td>';
      tr.onclick = () => select(s.id); tb.appendChild(tr);
    }
    note.textContent = r.formula + ' · ' + r.nHobserved + ' of ' + r.nHtotal + ' H shown' + (r.nHobserved < r.nHtotal ? ' (rest are exchangeable O–H / N–H)' : '') + '. Empirical estimate: typically within 0.3 ppm.' + (r.h.signals.some(q => q.kind === 'vinyl') ? ' Alkene cis/trans (and J) follow the drawn geometry.' : '') + (r.notes.length ? ' ' + r.notes.join(' ') : '');
  } else {
    th.innerHTML = '<tr><th class="n">δ (ppm)</th><th class="n">C</th><th>Type</th><th>Atoms</th></tr>';
    for (const s of r.c.signals) {
      if (S.nuc === 'DEPT' && s.type === 'C') continue;
      const tr = document.createElement('tr'); tr.className = 'row'; tr.dataset.id = s.id;
      tr.innerHTML = '<td class="n">' + s.shift.toFixed(1) + '</td><td class="n">' + s.n + '</td><td>' + s.type + (S.nuc === 'DEPT' ? (s.type === 'CH2' ? ' (down)' : ' (up)') : '') + '</td><td class="n">' + s.atoms.map(a => a + 1).join(', ') + '</td>';
      tr.onclick = () => select(s.id); tb.appendChild(tr);
    }
    note.textContent = r.formula + ' · ' + r.c.signals.length + ' distinct carbon environments. Empirical estimate: typically within 5 ppm; peak heights are indicative only.' + (r.notes.length ? ' ' + r.notes.join(' ') : '');
  }
  markRow();
}
function markRow() { document.querySelectorAll('#tb tr').forEach(t => t.classList.toggle('sel', S.sel != null && +t.dataset.id === S.sel)); const t = document.querySelector('#tb tr.sel'); if (t) t.scrollIntoView({ block: 'nearest' }); }
function select(id) { S.sel = id; markRow(); drawStructure(); draw(); }

/* ---------------- structure view ---------------- */
function drawStructure() {
  const svg = $('#str'); svg.innerHTML = '';
  if (!S.res) return; const G = S.res.G, A = G.atoms;
  let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9; A.forEach(a => { x0 = Math.min(x0, a.x); x1 = Math.max(x1, a.x); y0 = Math.min(y0, -a.y); y1 = Math.max(y1, -a.y); });
  const L = Math.max(1, ...G.bonds.map(b => Math.hypot(A[b.a].x - A[b.b].x, A[b.a].y - A[b.b].y))) * 0, bl = (() => { const d = G.bonds.map(b => Math.hypot(A[b.a].x - A[b.b].x, A[b.a].y - A[b.b].y)).sort((a, b) => a - b); return d.length ? d[Math.floor(d.length / 2)] : 1; })();
  const m = bl * 0.9, vx = x0 - m, vy = y0 - m, vw = Math.max(bl * 2, x1 - x0) + 2 * m, vh = Math.max(bl * 2, y1 - y0) + 2 * m;
  svg.setAttribute('viewBox', [vx, vy, vw, vh].join(' '));
  const fg = css('--fg'), accent = css('--accent'), sw = bl * 0.075, ns = 'http://www.w3.org/2000/svg';
  const el = (n, a, p) => { const e = document.createElementNS(ns, n); for (const k in a) e.setAttribute(k, a[k]); (p || svg).appendChild(e); return e; };
  const sigs = signalsNow(), sel = sigs.find(s => s.id === S.sel), hi = new Set(sel ? sel.atoms : []);
  const colors = { O: '#E5534B', N: '#5B8DEF', S: '#C9A227', F: '#4DB87A', Cl: '#4DB87A', Br: '#B5603A', I: '#9B59B6', P: '#E08A2E' };
  const lab = i => A[i].el !== 'C' || A[i].nb.length === 0;
  for (const b of G.bonds) {
    const a = A[b.a], c = A[b.b]; let ax = a.x, ay = -a.y, cx = c.x, cy = -c.y; const dx = cx - ax, dy = cy - ay, len = Math.hypot(dx, dy) || 1, ux = dx / len, uy = dy / len, px = -uy, py = ux;
    const sh = (i, s) => lab(i) ? bl * 0.22 * s : 0; ax += ux * sh(b.a, 1); ay += uy * sh(b.a, 1); cx -= ux * sh(b.b, 1); cy -= uy * sh(b.b, 1);
    const offs = b.o === 1 ? [0] : b.o === 2 ? [-0.5, 0.5] : [-1, 0, 1], g = bl * 0.09;
    for (const o of offs) el('line', { x1: ax + px * o * g * 1.6, y1: ay + py * o * g * 1.6, x2: cx + px * o * g * 1.6, y2: cy + py * o * g * 1.6, stroke: fg, 'stroke-width': sw, 'stroke-linecap': 'round' });
  }
  A.forEach((a, i) => {
    const x = a.x, y = -a.y;
    if (hi.has(i)) el('circle', { cx: x, cy: y, r: bl * 0.34, fill: accent, opacity: 0.35 });
    if (lab(i)) { const t = el('text', { x, y: y + bl * 0.15, 'text-anchor': 'middle', 'font-size': bl * 0.44, 'font-family': 'IBM Plex Sans, sans-serif', 'font-weight': 600, fill: colors[a.el] || fg }); t.textContent = a.el + (a.nH ? 'H' + (a.nH > 1 ? a.nH : '') : '') + (a.chg ? (a.chg > 0 ? '+' : '–') : ''); }
    if (S.nums) { const t = el('text', { x: x + bl * 0.2, y: y - bl * 0.18, 'font-size': bl * 0.27, fill: accent, 'font-family': 'IBM Plex Mono, monospace' }); t.textContent = i + 1; }
    const hit = el('circle', { cx: x, cy: y, r: bl * 0.34, fill: 'transparent', style: 'cursor:pointer' });
    hit.addEventListener('click', () => { const s = sigs.find(q => q.atoms.includes(i)); select(s ? s.id : null); });
  });
}

/* ---------------- controls ---------------- */
document.querySelectorAll('#nuc button').forEach(b => b.onclick = () => {
  document.querySelectorAll('#nuc button').forEach(x => x.classList.toggle('on', x === b)); S.nuc = b.dataset.v; S.sel = null; S.exps = []; S.armed = false; syncExpUI();
  $('#exwrap').hidden = S.nuc !== '1H'; $('#intwrap').hidden = S.nuc !== '1H'; buildSpectrum(true); renderTable(); drawStructure(); draw();
});
$('#mhz').onchange = e => { S.mhz = +e.target.value; buildSpectrum(false); draw(); };
$('#lw').oninput = e => { S.lw = +e.target.value; $('#lwv').textContent = S.lw.toFixed(1) + ' Hz'; buildSpectrum(false); draw(); };
$('#solv').onchange = e => { S.solv = e.target.value; if (S.last) update(); };
$('#ex').onchange = e => { S.ex = e.target.checked; buildSpectrum(true); renderTable(); draw(); };
$('#nums').onchange = e => { S.nums = e.target.checked; drawStructure(); };
$('#b-reset').onclick = () => { if (S.res) { S.view = autoView(); S.vz = 1; draw(); } };
$('#b-exp').onclick = () => { if (!S.y) return toast('Draw a structure first.'); S.armed = !S.armed; syncExpUI(); if (S.armed) toast('Drag across the spectrum to choose the region to expand.'); };
$('#b-expclr').onclick = () => { S.exps = []; syncExpUI(); draw(); };
$('#icurve').onchange = e => { S.opt.icurve = e.target.checked; draw(); };
$('#ibrace').onchange = e => { S.opt.brace = e.target.checked; draw(); };
$('#b-theme').onclick = () => window.CPTheme.toggle();

function report() {
  if (!S.res) return '';
  const sh = S.solv === 'D2O' || S.solv === 'CD3OD';
  if (S.nuc === '1H') {
    const parts = S.res.h.signals.filter(s => !s.exch || (S.ex && !sh)).map(s => {
      const j = s.Jlist && s.Jlist.length ? ', J = ' + s.Jlist.map(g => g.J.toFixed(1)).join(', ') + ' Hz' : '';
      return s.shift.toFixed(2) + ' (' + s.mult + j + ', ' + s.n + 'H)';
    });
    return '1H NMR (' + S.mhz + ' MHz, ' + S.solv + ') δ ' + parts.join(', ') + '.';
  }
  const sigs = S.res.c.signals.filter(s => S.nuc !== 'DEPT' || s.type !== 'C');
  return (S.nuc === 'DEPT' ? 'DEPT-135 ' : '13C NMR ') + '(' + (S.mhz * 0.2515).toFixed(0) + ' MHz, ' + S.solv + ') δ ' + sigs.map(s => s.shift.toFixed(1)).join(', ') + '.';
}
$('#b-copy').onclick = async () => { const t = report(); if (!t) return toast('Nothing to copy yet.'); try { await navigator.clipboard.writeText(t); toast('Copied.'); } catch (e) { const ta = document.createElement('textarea'); ta.value = t; document.body.appendChild(ta); ta.select(); try { document.execCommand('copy'); toast('Copied.'); } catch (x) { toast('Copy was blocked by the browser.'); } ta.remove(); } };
function save(name, blob) { const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }
$('#b-csv').onclick = () => {
  if (!S.res) return toast('Draw a structure first.');
  let t; if (S.nuc === '1H') t = 'shift_ppm,hydrogens,multiplicity,J_Hz,type,atoms\n' + S.res.h.signals.map(s => [s.shift.toFixed(3), s.n, s.mult, '"' + Jtext(s) + '"', '"' + envName(s) + '"', '"' + s.atoms.map(a => a + 1).join(' ') + '"'].join(',')).join('\n');
  else t = 'shift_ppm,carbons,type,atoms\n' + S.res.c.signals.map(s => [s.shift.toFixed(2), s.n, s.type, '"' + s.atoms.map(a => a + 1).join(' ') + '"'].join(',')).join('\n');
  save('compchemx-nmr-' + S.nuc + '.csv', new Blob([t + '\n'], { type: 'text/csv' }));
};
$('#b-png').onclick = () => {
  if (!S.y) return toast('Draw a structure first.');
  const bg = css('--stage'), c2 = document.createElement('canvas'); c2.width = cv.width; c2.height = cv.height; const g = c2.getContext('2d'); g.fillStyle = bg; g.fillRect(0, 0, c2.width, c2.height); g.drawImage(cv, 0, 0);
  c2.toBlob(b => save('compchemx-nmr-' + S.nuc + '.png', b), 'image/png');
};

/* ---------------- resizable panes ---------------- */
const LS = { get(k) { try { return localStorage.getItem('cp-nmr-' + k); } catch (e) { return null; } }, set(k, v) { try { localStorage.setItem('cp-nmr-' + k, v); } catch (e) {} } };
const main = $('#main'), right = $('#right');
const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
function setLW(px) { px = clamp(px, 300, Math.min(main.clientWidth - 420, 1100)); main.style.setProperty('--lw', px + 'px'); LS.set('lw', px); fitDraw(); draw(); }
function setLH(px) { px = clamp(px, 110, Math.max(150, right.clientHeight - $('#ctl').offsetHeight - 10 - 10 - 170)); right.style.setProperty('--lh', px + 'px'); LS.set('lh', px); draw(); }
function splitter(el, axis) {
  let start = 0, base = 0;
  el.addEventListener('pointerdown', e => { start = axis === 'v' ? e.clientX : e.clientY; base = axis === 'v' ? $('#left').offsetWidth : $('#low').offsetHeight; el.setPointerCapture(e.pointerId); el.classList.add('drag'); document.body.classList.add('dragging'); e.preventDefault(); });
  el.addEventListener('pointermove', e => { if (!el.classList.contains('drag')) return; if (axis === 'v') setLW(base + e.clientX - start); else setLH(base - (e.clientY - start)); });
  const end = () => { el.classList.remove('drag'); document.body.classList.remove('dragging'); };
  el.addEventListener('pointerup', end); el.addEventListener('pointercancel', end);
  el.addEventListener('dblclick', () => { if (axis === 'v') { main.style.removeProperty('--lw'); LS.set('lw', ''); fitDraw(); } else { right.style.removeProperty('--lh'); LS.set('lh', ''); } draw(); });
  el.addEventListener('keydown', e => { const k = e.key, d = e.shiftKey ? 60 : 20;
    if (axis === 'v' && (k === 'ArrowLeft' || k === 'ArrowRight')) { setLW($('#left').offsetWidth + (k === 'ArrowRight' ? d : -d)); e.preventDefault(); }
    if (axis === 'h' && (k === 'ArrowUp' || k === 'ArrowDown')) { setLH($('#low').offsetHeight + (k === 'ArrowUp' ? d : -d)); e.preventDefault(); } });
}
splitter($('#sp-v'), 'v'); splitter($('#sp-h'), 'h');
{ const w = +LS.get('lw'), h = +LS.get('lh'); if (w) main.style.setProperty('--lw', clamp(w, 300, 1100) + 'px'); if (h) right.style.setProperty('--lh', clamp(h, 110, 900) + 'px'); }
window.addEventListener('resize', () => { fitDraw(); });
$('#exwrap').hidden = false;
window.__nmr = { S, NMR, update, select };
loadDraw();
})();
