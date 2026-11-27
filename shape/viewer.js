/* viewer.js: three.js viewer for CCX Shape. Ball-and-stick, electron-domain clouds, polyhedron, dipole, angle labels. */
(function (root) {
'use strict';
const COL = { H: 0xffffff, C: 0x909090, N: 0x3050f8, O: 0xff2a2a, F: 0x90e050, Cl: 0x1ff01f, Br: 0xa62929, I: 0x940094, S: 0xffe030, P: 0xff8000, B: 0xffb5b5, Xe: 0x429eb0, Si: 0xf0c8a0, Al: 0xbfa6a6, Se: 0xffa100, Be: 0xc2ff00, Ge: 0x668f8f, As: 0xbd80e3, Sn: 0x668080, Sb: 0x9e63b5, Te: 0xd47a00, Kr: 0x5cb8d1, Ar: 0x80d1e3, Ne: 0xb3e3f5, He: 0xd9ffff };
const RAD = { H: .26, C: .38, N: .36, O: .35, F: .32, Cl: .5, Br: .55, I: .62, S: .5, P: .5, B: .4, Xe: .62 };

function create(host, THREE) {
  const T = THREE, scene = new T.Scene(), cam = new T.PerspectiveCamera(40, 1, .1, 200);
  const ren = new T.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
  ren.setPixelRatio(Math.min(2, window.devicePixelRatio || 1)); host.appendChild(ren.domElement); ren.domElement.style.display = 'block';
  scene.add(new T.AmbientLight(0xffffff, .65)); const dl = new T.DirectionalLight(0xffffff, .8); dl.position.set(3, 5, 6); scene.add(dl);
  const grp = new T.Group(); scene.add(grp);
  const Q0 = () => new T.Quaternion().setFromEuler(new T.Euler(.35, -.5, 0));
  grp.quaternion.copy(Q0());
  const turn = (ax, ay) => { const q = new T.Quaternion().setFromEuler(new T.Euler(ax, ay, 0, 'YXZ')); grp.quaternion.premultiply(q).normalize(); };
  let model = null, an = null, sel = -1, spin = false, dist = 9, radius = 3;
  const opt = { lp: true, bond: false, pi: false, en: false, poly: false, dip: true, bdip: true, ang: false, labels: true, style: 'ball' };
  const items = [], atomMeshes = [], v = (a) => new T.Vector3(a[0], a[1], a[2]);
  const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim() || '#888';
  function clear() { while (grp.children.length) { const o = grp.children.pop(); o.traverse(x => { if (x.geometry) x.geometry.dispose(); if (x.material) { if (x.material.map) x.material.map.dispose(); x.material.dispose(); } }); } atomMeshes.length = 0; }
  function sprite(text, col, heavy) {
    const c = document.createElement('canvas'), x = c.getContext('2d'); c.width = 256; c.height = 96; x.font = '600 44px IBM Plex Sans, system-ui, sans-serif'; x.textAlign = 'center'; x.textBaseline = 'middle';
    x.lineWidth = heavy ? 18 : 9; x.lineJoin = 'round'; x.strokeStyle = heavy ? '#15130f' : (css('--stage') || '#000'); x.strokeText(text, 128, 48); x.fillStyle = col; x.fillText(text, 128, 48);
    const m = new T.SpriteMaterial({ map: new T.CanvasTexture(c), depthTest: false, transparent: true }), s = new T.Sprite(m); s.scale.set(1.3, .49, 1); s.renderOrder = 10; return s;
  }
  function cyl(a, b, r, mat) { const d = b.clone().sub(a), L = d.length(), m = new T.Mesh(new T.CylinderGeometry(r, r, L, 14, 1, true), mat); m.position.copy(a.clone().add(b).multiplyScalar(.5)); m.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), d.normalize()); return m; }
  function lobe(pos, dir, len, wid, color, op) {   /* teardrop: lathe around +y, then aimed along dir */
    const pts = []; for (let k = 0; k <= 16; k++) { const t = k / 16, r = wid * Math.sin(Math.PI * t) * (1 - .35 * t); pts.push(new T.Vector2(r, len * t)); }
    const g = new T.LatheGeometry(pts, 24), m = new T.Mesh(g, new T.MeshPhongMaterial({ color, transparent: true, opacity: op, depthWrite: false, side: T.DoubleSide, shininess: 30 }));
    m.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), v(dir).normalize()); m.position.copy(pos); return m;
  }
  function build() {
    clear(); if (!model) return;
    const A = model.atoms, P = A.map(a => v([a.x, a.y, a.z])), bm = new T.MeshPhongMaterial({ color: css('--muted') === '#888' ? 0xaaaaaa : 0xb8b2a6, shininess: 40 });
    const blue = 0x4aa3ff, orange = 0xff9a4d, ball = opt.style === 'ball';
    A.forEach((a, i) => {
      const r = (RAD[a.el] || .45) * (ball ? 1 : .55) * (opt.style === 'wire' ? .3 : 1), m = new T.Mesh(new T.SphereGeometry(r, 24, 18), new T.MeshPhongMaterial({ color: COL[a.el] || 0xcc66cc, shininess: 60, emissive: i === sel ? 0x553322 : 0 }));
      m.position.copy(P[i]); m.userData.atom = i; grp.add(m); atomMeshes.push(m);
      const showEl = opt.labels && a.el !== 'H', showEN = opt.en && a.en && a.enShown !== false;
      if (showEl || showEN) { const chg = a.q ? (Math.abs(a.q) > 1 ? Math.abs(a.q) : '') + (a.q > 0 ? '+' : '−') : '', txt = (showEl ? a.el + chg : '') + (showEl && showEN ? ' ' : '') + (showEN ? a.en.toFixed(1) : '');
        const s = sprite(txt, showEN && !showEl ? '#ffd24d' : '#' + (COL[a.el] ? new T.Color(COL[a.el]).lerp(new T.Color(css('--fg')), .5).getHexString() : 'fff'), a.el === 'H'); s.position.copy(P[i]); s.scale.multiplyScalar(.7 * (txt.length > 4 ? 1.25 : 1)); grp.add(s); }
    });
    model.bonds.forEach(b => {
      const n = Math.max(1, Math.round(b.order)), a = P[b.a], c = P[b.b], d = c.clone().sub(a).normalize();
      let off = new T.Vector3(0, 1, 0).cross(d); if (off.lengthSq() < 1e-4) off = new T.Vector3(1, 0, 0); off.normalize();
      const rr = opt.style === 'wire' ? .03 : .075, frac = b.order % 1 > .1 && b.order % 1 < .9;
      for (let k = 0; k < n; k++) { const sh = (k - (n - 1) / 2) * .2; grp.add(cyl(a.clone().addScaledVector(off, sh), c.clone().addScaledVector(off, sh), rr, bm)); }
      if (frac) { const dm = new T.MeshPhongMaterial({ color: 0xb8b2a6, transparent: true, opacity: .4 }); grp.add(cyl(a.clone().addScaledVector(off, .2), c.clone().addScaledVector(off, .2), rr * .8, dm)); }
    });
    if (opt.lp) A.forEach((a, i) => {
      if (!(an && an[i] && an[i].center)) return;
      (a.lp || []).forEach(d => grp.add(lobe(P[i].clone().addScaledVector(v(d).normalize(), (RAD[a.el] || .4) * .55), d, .95, .42, orange, .38)));
      (a.rad || []).forEach(d => grp.add(lobe(P[i].clone().addScaledVector(v(d).normalize(), (RAD[a.el] || .4) * .55), d, .8, .3, 0xffdd55, .45)));
    });
    if (opt.bond && an) an.forEach(at => { if (!at.center) return; const i = at.i; model.bonds.forEach(b => { if (b.a !== i && b.b !== i) return; const o = b.a === i ? b.b : b.a, d = P[o].clone().sub(P[i]), L = d.length(); grp.add(lobe(P[i].clone().addScaledVector(d.clone().normalize(), .3), [d.x, d.y, d.z], L * .62, .3, blue, .22)); }); });
    if (opt.pi) model.bonds.forEach(b => {                      /* pi clouds above and below the bond plane */
      if (b.order < 1.45) return; const a = P[b.a], c = P[b.b], mid = a.clone().add(c).multiplyScalar(.5), d = c.clone().sub(a), L = d.length(); d.normalize();
      const nb = [...model.bonds.filter(x => x !== b && (x.a === b.a || x.b === b.a || x.a === b.b || x.b === b.b))]; let nrm = null;
      for (const x of nb) { const o = (x.a === b.a || x.a === b.b) ? x.b : x.a, w = P[o].clone().sub(P[x.a === o ? x.b : x.a]); const n2 = d.clone().cross(w); if (n2.lengthSq() > 1e-3) { nrm = n2.normalize(); break; } }
      if (!nrm) { nrm = new T.Vector3(0, 1, 0).cross(d); if (nrm.lengthSq() < 1e-4) nrm = new T.Vector3(1, 0, 0).cross(d); nrm.normalize(); }
      const planeN = nrm.clone().cross(d).normalize();           /* pi lobes sit perpendicular to the sigma plane */
      const lat = nrm.clone().cross(d).normalize(); void lat;
      [1, -1].forEach(sg => { const m = new T.Mesh(new T.SphereGeometry(1, 20, 14), new T.MeshPhongMaterial({ color: blue, transparent: true, opacity: .3, depthWrite: false })); m.scale.set(.32, L * .55, .32); m.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), d); m.position.copy(mid).addScaledVector(planeN, sg * .42); grp.add(m); });
      void planeN;
    });
    if (opt.poly && sel >= 0 && an && an[sel] && an[sel].SN >= 2) {
      const nb = []; model.bonds.forEach(b => { if (b.a === sel) nb.push(P[b.b].clone().sub(P[sel])); else if (b.b === sel) nb.push(P[b.a].clone().sub(P[sel])); });
      const pts = nb.map(x => x.clone().normalize().multiplyScalar(1.5)).concat((A[sel].lp || []).map(d => v(d).normalize().multiplyScalar(1.5)));
      const ln = new T.LineBasicMaterial({ color: new T.Color(css('--accent')), transparent: true, opacity: .75 }); const gm = [];
      for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++) { const d = pts[i].distanceTo(pts[j]); let mn = 1e9; for (let k = 0; k < pts.length; k++) if (k !== i) mn = Math.min(mn, pts[i].distanceTo(pts[k])); if (d < mn * 1.12) gm.push(pts[i], pts[j]); }
      const l = new T.LineSegments(new T.BufferGeometry().setFromPoints(gm), ln); l.position.copy(P[sel]); grp.add(l);
    }
    if (opt.ang && sel >= 0 && an && an[sel]) {
      an[sel].angles.slice(0, 12).forEach(g => { const u = P[g.i].clone().sub(P[sel]).normalize(), w = P[g.k].clone().sub(P[sel]).normalize(), m = u.clone().add(w); if (m.length() < 1e-3) return; m.normalize();
        const s = sprite(g.deg.toFixed(1) + '°', css('--accent')); s.position.copy(P[sel]).addScaledVector(m, 1.05); grp.add(s);
        const arc = []; for (let t = 0; t <= 12; t++) { const q = u.clone().lerp(w, t / 12).normalize().multiplyScalar(.7); arc.push(q.add(P[sel])); } grp.add(new T.Line(new T.BufferGeometry().setFromPoints(arc), new T.LineBasicMaterial({ color: new T.Color(css('--accent')) }))); });
    }
    const dp = an && an.dipole;
    const arrow = (a, b, rad, color, head) => { const mat = new T.MeshPhongMaterial({ color }), d = b.clone().sub(a), L = d.length(); if (L < 1e-3) return; d.normalize();
      const tip = b.clone().addScaledVector(d, -head * 1.9); grp.add(cyl(a, tip, rad, mat)); const h = new T.Mesh(new T.ConeGeometry(head, head * 2.2, 16), mat); h.position.copy(b.clone().addScaledVector(d, -head * 1.1)); h.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), d); grp.add(h); };
    if (opt.bdip && dp) dp.bonds.forEach(bd => {          /* one pink arrow per polar bond, beside the bond, pointing to the delta-minus end */
      if (bd.d < .5) return; const a = P[bd.a], c = P[bd.b], d = c.clone().sub(a), L = d.length(); d.normalize();
      let off = new T.Vector3(0, 1, 0).cross(d); if (off.lengthSq() < 1e-4) off = new T.Vector3(1, 0, 0).cross(d); off.normalize();
      const ra = (RAD[A[bd.a].el] || .4), rc = (RAD[A[bd.b].el] || .4), s0 = Math.min(L * .3, ra + .12), s1 = Math.min(L * .3, rc + .12), o = off.multiplyScalar(.42);
      arrow(a.clone().addScaledVector(d, s0).add(o), c.clone().addScaledVector(d, -s1).add(o), .035, 0xff4d9a, .1);
      const tail = new T.Mesh(new T.SphereGeometry(.07, 10, 8), new T.MeshPhongMaterial({ color: 0xff4d9a })); tail.position.copy(a.clone().addScaledVector(d, s0).add(o)); grp.add(tail);
      const dm = sprite('δ+', '#ff4d9a'), dn = sprite('δ−', '#ff4d9a'); dm.scale.multiplyScalar(.7); dn.scale.multiplyScalar(.7);
      dm.position.copy(a.clone().addScaledVector(d, Math.max(s0 - .1, .05)).addScaledVector(o, 1.5)); dn.position.copy(c.clone().addScaledVector(d, -Math.max(s1 - .1, .05)).addScaledVector(o, 1.5)); grp.add(dm); grp.add(dn); });
    if (opt.dip && dp && dp.polar) {                      /* net dipole in green, drawn clear of the atoms */
      const d = v(dp.vector).normalize(), c = new T.Vector3(); P.forEach(p => c.add(p)); c.multiplyScalar(1 / P.length);
      let perp = new T.Vector3(0, 1, 0).cross(d); if (perp.lengthSq() < 1e-4) perp = new T.Vector3(1, 0, 0).cross(d); perp.normalize();
      let best = perp, bx = 1e9; for (let k = 0; k < 12; k++) { const q = perp.clone().applyAxisAngle(d, k * Math.PI / 6); let e = 0; P.forEach(p => e = Math.max(e, p.clone().sub(c).dot(q))); if (e < bx) { bx = e; best = q; } }
      let lo = 1e9, hi = -1e9; P.forEach(p => { const t = p.clone().sub(c).dot(d); lo = Math.min(lo, t); hi = Math.max(hi, t); });
      const mid = c.clone().addScaledVector(d, (lo + hi) / 2).addScaledVector(best, bx + .95), L = Math.max(2.2, hi - lo + .9);
      arrow(mid.clone().addScaledVector(d, -L / 2), mid.clone().addScaledVector(d, L / 2), .07, 0x38d17a, .2);
      const s = sprite('net μ', '#38d17a'); s.position.copy(mid).addScaledVector(best, .5); grp.add(s);
    }
    radius = P.reduce((m, p) => Math.max(m, p.length()), 1) + 1.6; dist = Math.max(6, radius * 1.9);
  }
  function size() { const w = host.clientWidth || 300, h = host.clientHeight || 300; ren.setSize(w, h, false); ren.domElement.style.width = w + 'px'; ren.domElement.style.height = h + 'px'; cam.aspect = w / h; cam.updateProjectionMatrix(); }
  function frame() { cam.position.set(0, 0, dist); cam.lookAt(0, 0, 0); ren.render(scene, cam); }
  let raf = 0, last = 0, alive = true;
  function loop(t) { if (!alive) return; raf = requestAnimationFrame(loop); if (spin) turn(0, (t - last) * .0006); last = t; frame(); }
  raf = requestAnimationFrame(loop);
  /* interaction: drag rotate, wheel zoom, click select */
  let down = null, moved = 0; const el = ren.domElement; el.style.touchAction = 'none';
  el.addEventListener('pointerdown', e => { down = { x: e.clientX, y: e.clientY }; moved = 0; el.setPointerCapture(e.pointerId); });
  el.addEventListener('pointermove', e => { if (!down) return; const dx = e.clientX - down.x, dy = e.clientY - down.y; moved += Math.abs(dx) + Math.abs(dy); turn(dy * .01, dx * .01); down = { x: e.clientX, y: e.clientY }; });
  el.addEventListener('pointerup', e => { if (moved < 5) pick(e); down = null; });
  el.addEventListener('wheel', e => { e.preventDefault(); dist = Math.max(3, Math.min(60, dist * (e.deltaY > 0 ? 1.08 : .92))); }, { passive: false });
  const ray = new T.Raycaster();
  function pick(e) { const r = el.getBoundingClientRect(); ray.setFromCamera(new T.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1), cam); grp.updateMatrixWorld(true); const h = ray.intersectObjects(atomMeshes, false)[0]; api.onSelect && api.onSelect(h ? h.object.userData.atom : -1); }
  const ro = typeof ResizeObserver !== 'undefined' ? new ResizeObserver(size) : null; if (ro) ro.observe(host); size();
  const api = {
    load(m, analysis) { model = m; an = analysis ? analysis.atoms : null; if (an) an.dipole = analysis.polarity; sel = -1; build(); },
    select(i) { sel = i; build(); }, set(k, val) { opt[k] = val; build(); }, setSpin(b) { spin = b; }, reset() { grp.quaternion.copy(Q0()); dist = Math.max(6, radius * 1.9); },
    rebuild: build, resize: size, opt, png() { frame(); return ren.domElement.toDataURL('image/png'); },
    dispose() { alive = false; cancelAnimationFrame(raf); if (ro) ro.disconnect(); clear(); ren.dispose(); el.remove(); },
  };
  return api;
}
root.ShapeViewer = { create };
})(typeof self !== 'undefined' ? self : this);
