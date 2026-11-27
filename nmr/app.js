/* CompChemX NMR page logic: embeds CompChemX Draw, predicts with nmr.js, draws the spectrum, table and structure. */
(function () {
'use strict';
const $ = s => document.querySelector(s);
const S = { nuc: '1H', mhz: 400, lw: 1, solv: 'CDCl3', ex: true, nums: false, res: null, err: '', lines: [], y: null, grid: null, view: null, sel: null, last: null, drag: null, hover: null };
const dark = () => { const t = document.documentElement.dataset.theme; return t ? t === 'dark' : matchMedia('(prefers-color-scheme:dark)').matches; };
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const toast = m => { const t = $('#toast'); t.textContent = m; t.hidden = false; clearTimeout(toast.t); toast.t = setTimeout(() => t.hidden = true, 2200); };

/* ---------------- Draw bridge (same contract as CompChemX Dock: global molfile()) ---------------- */
const BL = { frame: $('#bl') };
const blWin = () => { try { return BL.frame.contentWindow; } catch (e) { return null; } };
function molfile() { const w = blWin(); try { if (w && typeof w.molfile === 'function') return w.molfile(); } catch (e) {} return null; }
function syncTheme() { const w = blWin(); try { if (w && w.document) w.document.documentElement.dataset.theme = dark() ? 'dark' : 'light'; } catch (e) {} }
async function loadDraw() {
  BL.frame.addEventListener('load', () => { syncTheme(); poll(); });
  try { const r = await fetch('../draw/index.html'); if (!r.ok) throw new Error('HTTP ' + r.status); BL.frame.srcdoc = await r.text(); }
  catch (e) { $('#msg').textContent = 'The CompChemX Draw editor could not be loaded (' + e.message + ').'; }
}
function poll() { if (S.hold) return; const t = molfile(); if (t == null || t === S.last) return; S.last = t; update(); }
setInterval(poll, 600);
window.addEventListener('cp-theme', () => { syncTheme(); draw(); drawStructure(); });

/* ---------------- prediction ---------------- */
function update() {
  S.sel = null; S.res = null; S.err = '';
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
const PAD = { l: 14, r: 14, t: 26, b: 48 };
function size() { const b = $('#specbox'), d = Math.min(2, window.devicePixelRatio || 1); const w = b.clientWidth, h = b.clientHeight; if (cv.width !== w * d || cv.height !== h * d) { cv.width = w * d; cv.height = h * d; } ctx.setTransform(d, 0, 0, d, 0, 0); return [w, h]; }
const xOf = (ppm, w) => PAD.l + (S.view[1] - ppm) / (S.view[1] - S.view[0]) * (w - PAD.l - PAD.r);
const ppmOf = (px, w) => S.view[1] - (px - PAD.l) / (w - PAD.l - PAD.r) * (S.view[1] - S.view[0]);
function niceStep(range, target) { const raw = range / target, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p; return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p; }
function signalsNow() { return !S.res ? [] : S.nuc === '1H' ? S.res.h.signals : S.res.c.signals; }

function draw() {
  const [w, h] = size(); ctx.clearRect(0, 0, w, h);
  const fg = css('--fg'), muted = css('--muted'), line = css('--line'), accent = css('--accent'), trace = css('--trace');
  if (!S.y || !S.view) return;
  const base = h - PAD.b, top = PAD.t, ph = base - top;
  const [v0, v1] = S.view, dx = (S.grid.hi - S.grid.lo) / (S.grid.N - 1);
  /* y scale: tallest visible point fills 82% of the plot (DEPT keeps zero on a baseline lifted for negatives) */
  const i0 = Math.max(0, Math.floor((v0 - S.grid.lo) / dx)), i1 = Math.min(S.grid.N - 1, Math.ceil((v1 - S.grid.lo) / dx));
  let mx = 0, mn = 0; for (let i = i0; i <= i1; i++) { const v = S.y[i]; if (v > mx) mx = v; if (v < mn) mn = v; }
  mx = mx || 1; const hasNeg = mn < -mx * 0.02, zero = hasNeg ? base - ph * 0.45 : base, up = hasNeg ? ph * 0.40 : ph * 0.80, sc = up / Math.max(mx, -mn * (hasNeg ? 1 : 0), 1e-12);
  const scale = hasNeg ? up / Math.max(mx, -mn) : up / mx;
  /* selection band */
  const sigs = signalsNow(), sel = sigs.find(s => s.id === S.sel);
  if (sel) {
    let a = sel.shift, b = sel.shift, pad = S.nuc === '1H' ? 0.04 : 1.5;
    if (S.nuc === '1H') { sel.lines.forEach(l => { const p = sel.shift + l.d / S.mhz; a = Math.min(a, p); b = Math.max(b, p); }); pad = 0.05 + (sel.exch ? 0.3 : 0); }
    ctx.fillStyle = css('--sel'); const xa = xOf(b + pad, w), xb = xOf(a - pad, w); ctx.fillRect(xa, top - 4, xb - xa, base - top + 4);
  }
  /* axis */
  ctx.strokeStyle = line; ctx.fillStyle = muted; ctx.lineWidth = 1; ctx.font = '11px ' + css('--mono'); ctx.textAlign = 'center';
  ctx.beginPath(); ctx.moveTo(PAD.l, base + 0.5); ctx.lineTo(w - PAD.r, base + 0.5); ctx.stroke();
  const st = niceStep(v1 - v0, 9), minor = st / 5;
  for (let t = Math.ceil(v0 / minor) * minor; t <= v1 + 1e-9; t += minor) { const x = Math.round(xOf(t, w)) + 0.5, big = Math.abs(t / st - Math.round(t / st)) < 1e-6; ctx.beginPath(); ctx.moveTo(x, base); ctx.lineTo(x, base + (big ? 6 : 3)); ctx.stroke(); if (big) ctx.fillText((+t.toFixed(3)).toString(), x, base + 19); }
  ctx.textAlign = 'right'; ctx.fillText('δ (ppm)', w - PAD.r, 14);
  ctx.textAlign = 'left'; ctx.fillText((S.nuc === '1H' ? '¹H ' + S.mhz + ' MHz' : S.nuc === '13C' ? '¹³C ' + (S.mhz * 0.2515).toFixed(0) + ' MHz' : 'DEPT-135 ' + (S.mhz * 0.2515).toFixed(0) + ' MHz') + ' · ' + S.solv + ' · simulated', PAD.l, 14);
  /* trace */
  ctx.save(); ctx.beginPath(); ctx.rect(PAD.l, 0, w - PAD.l - PAD.r, h); ctx.clip();
  ctx.strokeStyle = trace; ctx.lineWidth = 1.3; ctx.beginPath();
  const step = Math.max(1, Math.floor((i1 - i0) / (w * 2)));
  for (let i = i0; i <= i1; i += step) { let v = S.y[i]; if (step > 1) { for (let k = 1; k < step && i + k <= i1; k++) { const q = S.y[i + k]; if (Math.abs(q) > Math.abs(v)) v = q; } } const x = xOf(S.grid.lo + i * dx, w), y = zero - v * scale; if (i === i0) ctx.moveTo(x, y); else ctx.lineTo(x, y); }
  ctx.stroke();
  if (sel) { /* redraw the selected signal in the accent colour */
    ctx.strokeStyle = accent; ctx.lineWidth = 1.8; ctx.beginPath(); let first = true;
    const a = sel.shift - (S.nuc === '1H' ? 0.06 + (sel.exch ? 0.35 : 0) + (sel.lines[sel.lines.length - 1].d - sel.lines[0].d) / 2 / S.mhz : 1.2), b = sel.shift + (S.nuc === '1H' ? 0.06 + (sel.exch ? 0.35 : 0) + (sel.lines[sel.lines.length - 1].d - sel.lines[0].d) / 2 / S.mhz : 1.2);
    for (let i = Math.max(i0, Math.floor((a - S.grid.lo) / dx)); i <= Math.min(i1, Math.ceil((b - S.grid.lo) / dx)); i++) { const x = xOf(S.grid.lo + i * dx, w), y = zero - S.y[i] * scale; if (first) { ctx.moveTo(x, y); first = false; } else ctx.lineTo(x, y); }
    ctx.stroke();
  }
  ctx.restore();
  /* labels and integrals */
  ctx.font = '11px ' + css('--mono'); ctx.textAlign = 'center';
  const placed = [];
  const order = sigs.slice().sort((p, q) => q.shift - p.shift);
  for (const s of order) {
    if (S.nuc === '1H' && s.exch && !S.ex) continue; if (S.nuc === '1H' && s.exch && (S.solv === 'D2O' || S.solv === 'CD3OD')) continue;
    if (S.nuc === 'DEPT' && s.type === 'C') continue;
    if (s.shift < v0 || s.shift > v1) continue;
    const x = xOf(s.shift, w); let peak = 0; const ci = Math.round((s.shift - S.grid.lo) / dx);
    for (let k = -Math.round(0.15 / dx); k <= Math.round(0.15 / dx); k++) { const v = S.y[ci + k]; if (v !== undefined && Math.abs(v) > Math.abs(peak)) peak = v; }
    const ty = peak >= 0 ? zero - peak * scale - 5 : zero - peak * scale + 13; let lx = x;
    if (placed.some(p => Math.abs(p - lx) < 28)) continue; placed.push(lx);
    ctx.fillStyle = s.id === S.sel ? accent : fg; ctx.fillText(s.shift.toFixed(S.nuc === '1H' ? 2 : 1), x, Math.max(top - 2, ty));
  }
  if (S.nuc === '1H') { /* integration steps: one number per signal under the axis */
    ctx.fillStyle = accent; ctx.textAlign = 'center';
    for (const s of sigs) { if (s.exch && (!S.ex || S.solv === 'D2O' || S.solv === 'CD3OD')) continue; if (s.shift < v0 || s.shift > v1) continue; const x = xOf(s.shift, w); ctx.fillText(s.n + 'H', x, base + 35); }
  }
  S.geom = { w, h, zero, scale, base, top };
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
cv.addEventListener('mousedown', e => { S.drag = { x0: e.offsetX, x1: e.offsetX, moved: false }; });
window.addEventListener('mouseup', e => {
  const d = S.drag; S.drag = null; if (!d) return;
  if (d.moved && Math.abs(d.x1 - d.x0) > 6) { const w = cv.clientWidth, a = ppmOf(d.x0, w), b = ppmOf(d.x1, w); S.view = [Math.min(a, b), Math.max(a, b)]; draw(); }
  else if (!d.moved) { const s = sigAt(d.x0); select(s ? s.id : null); }
  $('#tip').hidden = true;
});
cv.addEventListener('mousemove', e => {
  if (S.drag) { S.drag.x1 = e.offsetX; if (Math.abs(S.drag.x1 - S.drag.x0) > 3) S.drag.moved = true; draw(); if (S.drag.moved) { ctx.fillStyle = css('--accent-soft'); ctx.fillRect(Math.min(S.drag.x0, S.drag.x1), PAD.t, Math.abs(S.drag.x1 - S.drag.x0), cv.clientHeight - PAD.t - PAD.b); } return; }
  if (!S.y) return; const s = sigAt(e.offsetX), tip = $('#tip'), p = ppmOf(e.offsetX, cv.clientWidth);
  tip.hidden = false; tip.textContent = p.toFixed(S.nuc === '1H' ? 3 : 1) + ' ppm' + (s ? '  ·  ' + (S.nuc === '1H' ? s.n + 'H ' + s.mult : s.n + ' C ' + s.type) : '');
  tip.style.left = Math.min(cv.clientWidth - tip.offsetWidth - 6, e.offsetX + 12) + 'px'; tip.style.top = Math.max(4, e.offsetY - 28) + 'px'; cv.style.cursor = s ? 'pointer' : 'crosshair';
});
cv.addEventListener('mouseleave', () => { $('#tip').hidden = true; });
cv.addEventListener('dblclick', () => { S.view = autoView(); draw(); });
cv.addEventListener('wheel', e => {
  if (!S.y) return; e.preventDefault(); const w = cv.clientWidth, p = ppmOf(e.offsetX, w), f = e.deltaY < 0 ? 0.8 : 1.25, [a, b] = S.view, [lo, hi] = full();
  let na = p - (p - a) * f, nb = p + (b - p) * f; if (nb - na > hi - lo) { na = lo; nb = hi; } S.view = [Math.max(lo, na), Math.min(hi, nb)]; draw();
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
  document.querySelectorAll('#nuc button').forEach(x => x.classList.toggle('on', x === b)); S.nuc = b.dataset.v; S.sel = null;
  $('#exwrap').hidden = S.nuc !== '1H'; buildSpectrum(true); renderTable(); drawStructure(); draw();
});
$('#mhz').onchange = e => { S.mhz = +e.target.value; buildSpectrum(false); draw(); };
$('#lw').oninput = e => { S.lw = +e.target.value; $('#lwv').textContent = S.lw.toFixed(1) + ' Hz'; buildSpectrum(false); draw(); };
$('#solv').onchange = e => { S.solv = e.target.value; if (S.last) update(); };
$('#ex').onchange = e => { S.ex = e.target.checked; buildSpectrum(true); renderTable(); draw(); };
$('#nums').onchange = e => { S.nums = e.target.checked; drawStructure(); };
$('#b-reset').onclick = () => { if (S.res) { S.view = autoView(); draw(); } };
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

$('#exwrap').hidden = false;
window.__nmr = { S, NMR, update, select };
loadDraw();
})();
