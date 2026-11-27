/* Shared navigation: a way back to the home page, and (where a page has no theme button of its own) a theme switch. */
(function () {
  if (window.parent !== window) return; /* embedded (Draw inside Dock): no navigation chrome */
  var home = new URL('../', location.href).href;
  var css = document.createElement('style');
  css.textContent = '.cp-nav{display:flex;align-items:center;gap:6px;margin-left:auto}' +
    '.cp-nav a,.cp-nav button{font:500 12.5px/1 var(--font-ui,var(--ui,var(--sans,system-ui,sans-serif)));height:30px;padding:0 11px;display:inline-flex;align-items:center;gap:6px;text-decoration:none;cursor:pointer;' +
    'color:var(--ink,var(--fg,inherit));background:var(--panel,transparent);border:1px solid var(--line,#444);border-radius:6px}' +
    '.cp-nav a:hover,.cp-nav button:hover{border-color:var(--accent)}' +
    '.cp-ver{font:500 11.5px/1 var(--font-mono,var(--mono,ui-monospace,monospace));color:var(--muted);border:1px solid var(--line,#444);border-radius:10px;padding:4px 8px;white-space:nowrap;align-self:center;display:inline-flex}';
  document.head.appendChild(css);
  var brand = document.querySelector('.brand');
  if (brand && brand.tagName !== 'A') {
    var a = document.createElement('a'); a.className = brand.className; a.href = home; a.title = 'All tools'; a.innerHTML = brand.innerHTML;
    a.style.cssText = 'color:inherit;text-decoration:none'; brand.replaceWith(a);
  }
  /* Home button: first thing in the header, filled accent colour, with a house icon, so it can't be mistaken for the tool's own buttons */
  var hdr = document.getElementById('top') || document.querySelector('header'), hh = document.createElement('a');
  hh.href = home; hh.id = 'cp-home'; hh.title = 'Back to all CompChemX tools'; hh.setAttribute('aria-label', 'Home: all tools');
  hh.innerHTML = '<svg viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M3 11.5L12 4l9 7.5"/><path d="M5.5 10v10h13V10"/><path d="M10 20v-5.5h4V20"/></svg><span>Home</span>';
  hh.style.cssText = 'flex:none;display:inline-flex;align-items:center;gap:6px;height:30px;padding:0 13px 0 11px;margin-right:6px;border-radius:6px;background:var(--accent);color:var(--accent-fg,#fff);font:600 13px/1 var(--font-ui,var(--ui,var(--sans,system-ui,sans-serif)));text-decoration:none;white-space:nowrap;box-shadow:0 1px 0 rgba(0,0,0,.25)';
  hh.addEventListener('mouseenter', function () { hh.style.filter = 'brightness(1.1)'; }); hh.addEventListener('mouseleave', function () { hh.style.filter = ''; });
  if (hdr) hdr.insertBefore(hh, hdr.firstChild); else { hh.style.cssText += ';position:fixed;top:10px;left:10px;z-index:50'; document.body.appendChild(hh); }
  if (!document.getElementById('b-theme') && !document.getElementById('cp-theme-btn')) { /* tools without a theme button of their own */
    var host = hdr || document.body, n = document.createElement('div'); n.className = 'cp-nav';
    n.innerHTML = '<button type="button" id="cp-theme-btn" title="Switch dark / light">Theme</button>';
    host.appendChild(n); n.querySelector('button').addEventListener('click', function () { window.CPTheme.toggle(); });
  }
})();
