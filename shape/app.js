/* CCX Shape page logic: embeds CCX Draw, runs shape.js, draws the Lewis structure, the 3D view and the VSEPR table. */
(function () {
'use strict';
const $ = s => document.querySelector(s);
const S = { res: null, err: '', last: null, sel: -1, form: 0, hold: false, spin: false };
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const dark = () => { const t = document.documentElement.dataset.theme; return t ? t === 'dark' : matchMedia('(prefers-color-scheme:dark)').matches; };
const toast = m => { const t = $('#toast'); t.textContent = m; t.hidden = false; clearTimeout(toast.t); toast.t = setTimeout(() => t.hidden = true, 2200); };
const esc = s => String(s).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
const sub = n => String(n).replace(/\d/g, d => '₀₁₂₃₄₅₆₇₈₉'[d]);

/* ---------------- Draw bridge ---------------- */
const BL = { frame: $('#bl') };
const blWin = () => { try { return BL.frame.contentWindow; } catch (e) { return null; } };
function molfile() { const w = blWin(); try { if (w && typeof w.molfile === 'function') return w.molfile(); } catch (e) {} return null; }
function prepareDraw() {
  const w = blWin(); try {
    const d = w.document, st = d.createElement('style');
    st.textContent = '@media (max-width:880px){html,body{height:auto!important}#app{height:auto!important;min-height:0!important;grid-template-rows:auto auto 560px auto auto!important}}';
    d.head.appendChild(st);
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
  catch (e) { $('#vmsg').textContent = 'The CCX Draw editor could not be loaded (' + e.message + ').'; }
}
function poll() { if (S.hold) return; const t = molfile(); if (t == null || t === S.last) return; S.last = t; update(); }
setInterval(() => { poll(); fitDraw(); }, 600);

/* ---------------- 3D viewer ---------------- */
let V = null;
function viewer() {
  if (V) return V;
  try { V = ShapeViewer.create($('#v3'), THREE); V.onSelect = i => select(i); } catch (e) { $('#vmsg').textContent = '3D view unavailable: ' + e.message; $('#vmsg').hidden = false; }
  return V;
}

/* ---------------- analysis ---------------- */
function update() {
  const t = S.last; S.err = ''; S.res = null;
  try {
    if (!t || !/\n/.test(t)) throw new Error('Draw a structure to begin.');
    S.res = Shape.analyze(t, { chem: root.Chem });
    S.form = 0; const cs = S.res.atoms.filter(a => a.center); S.sel = cs.length ? cs.reduce((b, a) => a.SN > b.SN ? a : b).i : -1;
  } catch (e) { S.err = e.message || String(e); }
  render();
}
const root = window;
function select(i) { S.sel = S.sel === i ? -1 : i; if (V) V.select(S.sel); drawTable(); }

function render() {
  const msg = $('#vmsg'); msg.hidden = !S.err && !!S.res; msg.textContent = S.err || (S.res ? '' : 'Draw a molecule to begin.');
  if (!S.res) { $('#lew').innerHTML = ''; $('#tb').innerHTML = ''; $('#sum').innerHTML = ''; $('#notes').innerHTML = ''; $('#r-lab').textContent = 'Lewis structure'; const v = viewer(); if (v) v.load(null); return; }
  const v = viewer(); if (v) { v.load(S.res.model, S.res); v.select(S.sel); }
  drawLewis(); drawTable();
}

/* ---------------- Lewis diagram ---------------- */
function drawLewis() {
  const R = S.res, svg = $('#lew'), forms = R.resonance.forms, f = forms[S.form % forms.length], G = R.G, A = G.atoms;
  const lay = Shape.layout2D(G, f.state), P = lay.pos, L = lay.bl;
  let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9; P.forEach(p => { x0 = Math.min(x0, p.x); x1 = Math.max(x1, p.x); y0 = Math.min(y0, p.y); y1 = Math.max(y1, p.y); });
  lay.dots.forEach(p => { x0 = Math.min(x0, p.x); x1 = Math.max(x1, p.x); y0 = Math.min(y0, p.y); y1 = Math.max(y1, p.y); });
  const pad = L * 0.9, w = Math.max(x1 - x0 + 2 * pad, L * 2.6), h = Math.max(y1 - y0 + 2 * pad, L * 2.6), cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
  svg.setAttribute('viewBox', [cx - w / 2, cy - h / 2, w, h].join(' '));
  const fg = css('--fg'), muted = css('--muted'), acc = css('--accent'), bg = css('--paper'), fs = L * 0.34, parts = [];
  const gap = L * 0.2;
  G.bonds.forEach((b, k) => {
    const a = P[b.a], c = P[b.b], dx = c.x - a.x, dy = c.y - a.y, d = Math.hypot(dx, dy) || 1, ux = dx / d, uy = dy / d, nx = -uy, ny = ux, o = Math.round(f.state.bo[k]);
    const sa = (A[b.a].el !== 'C' || A[b.a].nb.length === 1 || true) ? gap : 0, sc = gap;
    for (let m = 0; m < o; m++) { const s = (m - (o - 1) / 2) * L * 0.12; parts.push(`<line x1="${a.x + ux * sa + nx * s}" y1="${a.y + uy * sa + ny * s}" x2="${c.x - ux * sc + nx * s}" y2="${c.y - uy * sc + ny * s}" stroke="${fg}" stroke-width="${L * 0.04}" stroke-linecap="round"/>`); }
  });
  A.forEach((a, i) => { const q = f.state.q[i], p = P[i];
    parts.push(`<text x="${p.x}" y="${p.y}" font-size="${fs * 1.15}" font-weight="600" text-anchor="middle" dominant-baseline="central" fill="${a.imp ? muted : fg}" font-family="IBM Plex Sans,system-ui,sans-serif">${a.el}</text>`);
    if (q) parts.push(`<g><circle cx="${p.x + fs * .9}" cy="${p.y - fs * .85}" r="${fs * .42}" fill="none" stroke="${q > 0 ? '#E6A33D' : acc}" stroke-width="${L * .025}"/><text x="${p.x + fs * .9}" y="${p.y - fs * .85}" font-size="${fs * .65}" text-anchor="middle" dominant-baseline="central" fill="${q > 0 ? '#E6A33D' : acc}">${q > 0 ? '+' : '−'}${Math.abs(q) > 1 ? Math.abs(q) : ''}</text></g>`); });
  lay.dots.forEach(d => parts.push(`<circle cx="${d.x}" cy="${d.y}" r="${L * .035}" fill="${acc}"/>`));
  void bg; svg.innerHTML = parts.join('');
  const n = forms.length; $('#r-lab').textContent = n > 1 ? `Resonance form ${S.form % n + 1} of ${n}${f.major ? ' · major' : ' · minor'}` : 'Lewis structure';
  $('#r-prev').disabled = $('#r-next').disabled = n < 2;
}

/* ---------------- table ---------------- */
function drawTable() {
  const R = S.res; if (!R) return; const A = R.G.atoms, tb = $('#tb');
  const lbl = i => A[i].el + (A[i].imp ? '' : (A[i].src + 1));
  tb.innerHTML = R.atoms.filter(a => a.center).map(a => {
    const angs = a.angles.length > 8 ? summarise(a.angles) : a.angles.map(g => `${lbl(g.i)}–${lbl(a.i)}–${lbl(g.k)} ${g.deg.toFixed(1)}°`).join(', ');
    return `<tr class="row${a.i === S.sel ? ' sel' : ''}" data-i="${a.i}"><td>${esc(lbl(a.i))}</td><td class="n">${a.axe.replace(/\d/g, d => sub(d)).replace('AX', 'AX').replace(/E(\d?)/, (m, d) => 'E' + sub(d || ''))}</td><td>${a.electronGeometry}</td><td><b>${a.shape}</b></td><td>${a.hybridisation || '–'}</td><td class="n" style="white-space:normal">${angs}</td></tr>`;
  }).join('') || '<tr><td colspan="6" style="color:var(--muted)">No central atoms: this is an atom or a diatomic with no bond angle.</td></tr>';
  tb.querySelectorAll('tr.row').forEach(tr => tr.addEventListener('click', () => select(+tr.dataset.i)));
  const Lw = R.lewis, pol = R.polarity;
  $('#sum').innerHTML = `<span><b>${Shape.formula(R.G).replace(/\d+/g, d => sub(d))}</b></span><span>Valence electrons <span class="chip">${Lw.totalValence}</span> (${Lw.pairs} pairs)</span><span>Bonding <span class="chip">${Lw.bondingElectrons}</span> · lone pair <span class="chip">${Lw.lonePairElectrons}</span></span><span>Net charge <span class="chip">${Lw.netCharge > 0 ? '+' : Lw.netCharge < 0 ? '−' : ''}${Math.abs(Lw.netCharge)}</span></span><span>${pol.polar ? '<b>Polar</b>' : pol.hasPolarBonds ? '<b>Nonpolar</b> (polar bonds cancel)' : '<b>Nonpolar</b>'} <span style="color:var(--muted)">· bond-dipole sum ${pol.magnitude.toFixed(2)} (ΔEN units)</span></span>`;
  const notes = Lw.notes.concat(R.warnings || []).map(n => `<div class="note${/exceeds|incomplete|more/.test(n) ? ' w' : ''}">${esc(n)}</div>`);
  notes.push(`<div class="note">3D model: ${esc(R.method)}. Angles come from a weighted electron-domain repulsion model in which lone pairs squeeze bond angles. VSEPR applies to main-group molecules, not transition-metal complexes.</div>`);
  $('#notes').innerHTML = notes.join('');
}
function summarise(angs) { const b = {}; angs.forEach(g => { const k = Math.round(g.deg); b[k] = (b[k] || 0) + 1; }); return Object.keys(b).sort((x, y) => x - y).map(k => `${k}° ×${b[k]}`).join(', '); }

/* ---------------- controls ---------------- */
$('#r-prev').addEventListener('click', () => { const n = S.res.resonance.forms.length; S.form = (S.form + n - 1) % n; drawLewis(); });
$('#r-next').addEventListener('click', () => { const n = S.res.resonance.forms.length; S.form = (S.form + 1) % n; drawLewis(); });
document.querySelectorAll('#ctl input[data-o]').forEach(c => c.addEventListener('change', () => { const v = viewer(); if (v) v.set(c.dataset.o, c.checked); }));
document.querySelectorAll('#sty button').forEach(b => b.addEventListener('click', () => { document.querySelectorAll('#sty button').forEach(x => x.classList.toggle('on', x === b)); const v = viewer(); if (v) v.set('style', b.dataset.v); }));
$('#b-spin').addEventListener('click', e => { S.spin = !S.spin; e.target.classList.toggle('on', S.spin); const v = viewer(); if (v) v.setSpin(S.spin); });
$('#b-reset').addEventListener('click', () => { const v = viewer(); if (v) v.reset(); });
$('#b-png').addEventListener('click', () => { const v = viewer(); if (!v || !S.res) return; const a = document.createElement('a'); a.href = v.png(); a.download = 'compchemx-shape.png'; a.click(); });
window.addEventListener('cp-theme', () => { syncTheme(); if (S.res) { drawLewis(); if (V) V.rebuild(); } });

/* ---------------- resizable panes ---------------- */
const LS = { get(k) { try { return localStorage.getItem('cp-shape-' + k); } catch (e) { return null; } }, set(k, v) { try { localStorage.setItem('cp-shape-' + k, v); } catch (e) {} } };
const main = $('#main'), right = $('#right');
const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
function setLW(px) { px = clamp(px, 300, Math.min(main.clientWidth - 420, 1100)); main.style.setProperty('--lw', px + 'px'); LS.set('lw', px); fitDraw(); if (V) V.resize(); }
function setLH(px) { px = clamp(px, 110, Math.max(150, right.clientHeight - $('#ctl').offsetHeight - 10 - 10 - 230)); right.style.setProperty('--lh', px + 'px'); LS.set('lh', px); if (V) V.resize(); }
function splitter(el, axis) {
  let start = 0, base = 0;
  el.addEventListener('pointerdown', e => { start = axis === 'v' ? e.clientX : e.clientY; base = axis === 'v' ? $('#left').offsetWidth : $('#anabox').offsetHeight; el.setPointerCapture(e.pointerId); el.classList.add('drag'); document.body.classList.add('dragging'); e.preventDefault(); });
  el.addEventListener('pointermove', e => { if (!el.classList.contains('drag')) return; if (axis === 'v') setLW(base + e.clientX - start); else setLH(base - (e.clientY - start)); });
  const end = () => { el.classList.remove('drag'); document.body.classList.remove('dragging'); };
  el.addEventListener('pointerup', end); el.addEventListener('pointercancel', end);
  el.addEventListener('dblclick', () => { if (axis === 'v') { main.style.removeProperty('--lw'); LS.set('lw', ''); fitDraw(); } else { right.style.removeProperty('--lh'); LS.set('lh', ''); } if (V) V.resize(); });
  el.addEventListener('keydown', e => { const k = e.key, d = e.shiftKey ? 60 : 20;
    if (axis === 'v' && (k === 'ArrowLeft' || k === 'ArrowRight')) { setLW($('#left').offsetWidth + (k === 'ArrowRight' ? d : -d)); e.preventDefault(); }
    if (axis === 'h' && (k === 'ArrowUp' || k === 'ArrowDown')) { setLH($('#anabox').offsetHeight + (k === 'ArrowUp' ? d : -d)); e.preventDefault(); } });
}
splitter($('#sp-v'), 'v'); splitter($('#sp-h'), 'h');
{ const w = +LS.get('lw'), h = +LS.get('lh'); if (w) main.style.setProperty('--lw', clamp(w, 300, 1100) + 'px'); if (h) right.style.setProperty('--lh', clamp(h, 110, 900) + 'px'); }
window.addEventListener('resize', () => { fitDraw(); if (V) V.resize(); });
window.__shape = { S, Shape, update, select };
loadDraw(); viewer();
})();
