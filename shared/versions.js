/* One place for every version number. Bump a tool's number when you change that tool, and "site" when you change the home page.
   Keys match the folder names. Everything is at 0.3.3 as of 27 Nov 2026. */
window.CP_VERSIONS = {
  site: '0.3.3',
  'draw': '0.3.3',
  'dock': '0.3.3',
  'nmr': '0.3.3'
};
/* Display text: a number such as "0.2.1" is shown as "v0.2.1". */
window.cpVer = function (v) { return /^\d/.test(v) ? 'v' + v : v; };
