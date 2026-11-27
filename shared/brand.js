/* One place for every product name. To rename the platform, change `platform` below (and nothing else).
   - The whole platform is `platform` (now "CompChemX").
   - Each tool is shown as "<platform> <suffix>", keyed by its folder name: draw/ is the 2D drawing tool, dock/ is the docking tool, nmr/ is the NMR simulator.
   Folder names, URLs and the version keys in versions.js do not change when you rename, so existing links keep working.
   Pages show these names through data-brand attributes (filled in when the page loads) and through window.cpName(). */
window.CP_BRAND = {
  platform: 'CompChemX',
  tools: { draw: 'Draw', dock: 'Dock', nmr: 'NMR', shape: 'Shape' }
};
/* cpName('site') -> platform name; cpName('draw') -> "CompChemX Draw"; cpName('dock') -> "CompChemX Dock" */
window.cpName = function (id) {
  var b = window.CP_BRAND;
  return !id || id === 'site' || id === 'platform' ? b.platform : (b.tools[id] ? b.platform + ' ' + b.tools[id] : b.platform);
};
/* <span data-brand="platform"></span>, <span data-brand="dock"></span>; optional data-brand-suffix=" Changelog" for titles. */
window.cpBrandApply = function (root) {
  (root || document).querySelectorAll('[data-brand]').forEach(function (e) { e.textContent = window.cpName(e.getAttribute('data-brand')) + (e.getAttribute('data-brand-suffix') || ''); });
};
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', function () { window.cpBrandApply(); }); else window.cpBrandApply();
