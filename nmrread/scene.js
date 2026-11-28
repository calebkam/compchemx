/* CCX NMRRead - drawing primitives and writers (SVG, PDF). Pure; works in browsers and node.
   A scene is a list of primitives in CSS pixels (96 per inch), y down:
     {t:'line', x1,y1,x2,y2, lw, c, dash}   {t:'path', pts:[[x,y],...], lw, c, dash}
     {t:'rect', x,y,w,h, fill, stroke, lw}  {t:'text', x,y, size, anchor:'start|middle|end', runs, c, rot, bold}
   Text runs: [{s:'7.26', k:'n'}], k = n (normal) | sup | sub | sym (Greek, drawn in the Symbol font in PDF). */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.NMRScene = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  /* "{^13}C NMR (400 MHz, CDCl{_3}) δ" -> runs */
  function rich(str) {
    const runs = []; let buf = '';
    const flush = () => { if (buf) { runs.push({ s: buf, k: 'n' }); buf = ''; } };
    for (let i = 0; i < str.length; i++) {
      const c = str[i];
      if (c === '{' && (str[i + 1] === '^' || str[i + 1] === '_')) {
        const j = str.indexOf('}', i); if (j > 0) { flush(); runs.push({ s: str.slice(i + 2, j), k: str[i + 1] === '^' ? 'sup' : 'sub' }); i = j; continue; }
      }
      if (c === 'δ') { flush(); runs.push({ s: 'δ', k: 'sym' }); continue; }
      buf += c;
    }
    flush();
    return runs;
  }
  const plain = runs => runs.map(r => r.s).join('');
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const f2 = n => (Math.round(n * 100) / 100).toString();

  /* width of a text primitive; mw(str, size, bold) measures plain text */
  function textWidth(p, mw) {
    let w = 0;
    for (const r of p.runs) {
      const sz = r.k === 'sup' || r.k === 'sub' ? p.size * 0.7 : p.size;
      w += r.k === 'sym' ? sz * 0.5 : mw(r.s, sz, p.bold);
    }
    return w;
  }

  /* ---------- SVG ---------- */
  function toSVG(prims, W, H, o) {
    o = o || {};
    const out = [];
    if (o.bg) out.push('<rect width="' + W + '" height="' + H + '" fill="' + o.bg + '"/>');
    const da = d => d && d.length ? ' stroke-dasharray="' + d.join(' ') + '"' : '';
    for (const p of prims) {
      if (p.t === 'line') out.push('<line x1="' + f2(p.x1) + '" y1="' + f2(p.y1) + '" x2="' + f2(p.x2) + '" y2="' + f2(p.y2) + '" stroke="' + p.c + '" stroke-width="' + f2(p.lw) + '"' + da(p.dash) + '/>');
      else if (p.t === 'path') out.push('<path d="' + p.pts.map((q, i) => (i ? 'L' : 'M') + f2(q[0]) + ' ' + f2(q[1])).join('') + '" fill="none" stroke="' + p.c + '" stroke-width="' + f2(p.lw) + '" stroke-linejoin="round" stroke-linecap="round"' + da(p.dash) + '/>');
      else if (p.t === 'rect') out.push('<rect x="' + f2(p.x) + '" y="' + f2(p.y) + '" width="' + f2(p.w) + '" height="' + f2(p.h) + '" fill="' + (p.fill || 'none') + '"' + (p.stroke ? ' stroke="' + p.stroke + '" stroke-width="' + f2(p.lw || 1) + '"' : '') + '/>');
      else if (p.t === 'text') {
        let prev = 0, tsp = '';
        for (const r of p.runs) {
          const rise = r.k === 'sup' ? p.size * 0.38 : r.k === 'sub' ? -p.size * 0.18 : 0, sz = r.k === 'sup' || r.k === 'sub' ? p.size * 0.7 : p.size;
          tsp += '<tspan' + (rise !== prev ? ' dy="' + f2(-(rise - prev)) + '"' : '') + (sz !== p.size ? ' font-size="' + f2(sz) + '"' : '') + '>' + esc(r.s) + '</tspan>';
          prev = rise;
        }
        const tr = p.rot ? ' transform="translate(' + f2(p.x) + ' ' + f2(p.y) + ') rotate(' + p.rot + ')"' : ' x="' + f2(p.x) + '" y="' + f2(p.y) + '"';
        out.push('<text' + tr + ' font-size="' + f2(p.size) + '"' + (p.bold ? ' font-weight="600"' : '') + ' text-anchor="' + (p.anchor || 'start') + '" fill="' + p.c + '">' + tsp + '</text>');
      }
    }
    return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' + W + ' ' + H + '" font-family="' + esc(o.font || 'Arial, Helvetica, sans-serif') + '"' + (o.mm ? ' width="' + o.mm[0] + 'mm" height="' + o.mm[1] + 'mm"' : '') + '>' + out.join('') + '</svg>';
  }

  /* ---------- PDF (vector, base-14 fonts) ---------- */
  const WIN = { '–': 0x96, '—': 0x97, '×': 0xD7, '°': 0xB0, '±': 0xB1, 'µ': 0xB5, '·': 0xB7, '‘': 0x91, '’': 0x92, '“': 0x93, '”': 0x94, '•': 0x95, '¹': 0xB9, '²': 0xB2, '³': 0xB3, '−': 0x2D };
  function pdfStr(s) {
    let o = '(';
    for (const ch of s) {
      let c = ch.charCodeAt(0);
      if (c > 255 || c < 32) c = WIN[ch] || 63; else if (WIN[ch]) c = WIN[ch];
      if (c === 40 || c === 41 || c === 92) o += '\\' + String.fromCharCode(c);
      else o += String.fromCharCode(c);
    }
    return o + ')';
  }
  function rgb(hex) {
    let h = String(hex || '#000000').replace('#', ''); if (h.length === 3) h = h.replace(/./g, '$&$&');
    const n = parseInt(h, 16) || 0;
    return [(n >> 16 & 255) / 255, (n >> 8 & 255) / 255, (n & 255) / 255].map(v => f2(v)).join(' ');
  }
  /* o: {mm:[w,h], mw(str,size,bold), serif:boolean, title, bg} -> Uint8Array */
  function toPDF(prims, W, H, o) {
    const k = 0.75, PW = f2(o.mm[0] * 72 / 25.4), PH = f2(o.mm[1] * 72 / 25.4), sx = o.mm[0] * 72 / 25.4 / W, sy = o.mm[1] * 72 / 25.4 / H;
    const X = x => f2(x * sx), Y = y => f2(o.mm[1] * 72 / 25.4 - y * sy);
    const fSans = o.serif ? '/F3' : '/F1', fBold = o.serif ? '/F5' : '/F2', c = [];
    if (o.bg) c.push(rgb(o.bg) + ' rg 0 0 ' + PW + ' ' + PH + ' re f');
    const stroke = p => { c.push(rgb(p.c) + ' RG ' + f2(p.lw * k) + ' w 1 J 1 j ' + (p.dash && p.dash.length ? '[' + p.dash.map(d => f2(d * k)).join(' ') + '] 0 d' : '[] 0 d')); };
    for (const p of prims) {
      if (p.t === 'line') { stroke(p); c.push(X(p.x1) + ' ' + Y(p.y1) + ' m ' + X(p.x2) + ' ' + Y(p.y2) + ' l S'); }
      else if (p.t === 'path') { stroke(p); c.push(p.pts.map((q, i) => X(q[0]) + ' ' + Y(q[1]) + (i ? ' l' : ' m')).join(' ') + ' S'); }
      else if (p.t === 'rect') {
        if (p.fill) c.push(rgb(p.fill) + ' rg ' + X(p.x) + ' ' + Y(p.y + p.h) + ' ' + f2(p.w * sx) + ' ' + f2(p.h * sy) + ' re f');
        if (p.stroke) { stroke({ c: p.stroke, lw: p.lw || 1 }); c.push(X(p.x) + ' ' + Y(p.y + p.h) + ' ' + f2(p.w * sx) + ' ' + f2(p.h * sy) + ' re S'); }
      } else if (p.t === 'text') {
        const w = textWidth(p, o.mw), ox = p.anchor === 'middle' ? -w / 2 : p.anchor === 'end' ? -w : 0;
        const th = -(p.rot || 0) * Math.PI / 180, cs = Math.cos(th), sn = Math.sin(th);
        c.push(rgb(p.c) + ' rg');
        let adv = ox;
        for (const r of p.runs) {
          const sz = r.k === 'sup' || r.k === 'sub' ? p.size * 0.7 : p.size, rise = r.k === 'sup' ? p.size * 0.38 : r.k === 'sub' ? -p.size * 0.18 : 0;
          /* position of this run in scene coordinates: along the baseline direction (cos, -sin in y-down space), raised perpendicular */
          const lx = adv, ly = -rise, px = p.x + lx * Math.cos((p.rot || 0) * Math.PI / 180) - ly * Math.sin((p.rot || 0) * Math.PI / 180), py = p.y + lx * Math.sin((p.rot || 0) * Math.PI / 180) + ly * Math.cos((p.rot || 0) * Math.PI / 180);
          const font = r.k === 'sym' ? '/F4' : (p.bold ? fBold : fSans), s = r.k === 'sym' ? 'd' : r.s;
          c.push('BT ' + font + ' ' + f2(sz * sx) + ' Tf ' + f2(cs) + ' ' + f2(sn) + ' ' + f2(-sn) + ' ' + f2(cs) + ' ' + X(px) + ' ' + Y(py) + ' Tm ' + (r.k === 'sym' ? '(d)' : pdfStr(s)) + ' Tj ET');
          adv += r.k === 'sym' ? sz * 0.5 : o.mw(r.s, sz, p.bold);
        }
      }
    }
    const stream = c.join('\n');
    const objs = [
      '<< /Type /Catalog /Pages 2 0 R >>',
      '<< /Type /Pages /Kids [3 0 R] /Count 1 >>',
      '<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ' + PW + ' ' + PH + '] /Contents 4 0 R /Resources << /Font << /F1 5 0 R /F2 6 0 R /F3 7 0 R /F4 8 0 R /F5 9 0 R >> >> >>',
      '<< /Length ' + stream.length + ' >>\nstream\n' + stream + '\nendstream',
      '<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>',
      '<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>',
      '<< /Type /Font /Subtype /Type1 /BaseFont /Times-Roman /Encoding /WinAnsiEncoding >>',
      '<< /Type /Font /Subtype /Type1 /BaseFont /Symbol >>',
      '<< /Type /Font /Subtype /Type1 /BaseFont /Times-Bold /Encoding /WinAnsiEncoding >>',
      '<< /Title ' + pdfStr(o.title || 'NMR spectrum') + ' /Creator (CCX NMRRead) >>'
    ];
    let pdf = '%PDF-1.4\n%\xE2\xE3\xCF\xD3\n'; const offs = [];
    objs.forEach((b, i) => { offs.push(pdf.length); pdf += (i + 1) + ' 0 obj\n' + b + '\nendobj\n'; });
    const xr = pdf.length;
    pdf += 'xref\n0 ' + (objs.length + 1) + '\n0000000000 65535 f \n' + offs.map(n => String(n).padStart(10, '0') + ' 00000 n \n').join('');
    pdf += 'trailer\n<< /Size ' + (objs.length + 1) + ' /Root 1 0 R /Info 10 0 R >>\nstartxref\n' + xr + '\n%%EOF\n';
    const u = new Uint8Array(pdf.length); for (let i = 0; i < pdf.length; i++) u[i] = pdf.charCodeAt(i) & 255;
    return u;
  }

  return { rich, plain, textWidth, toSVG, toPDF };
});
