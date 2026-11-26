/* One place for every version number. Bump a tool's number when you change that tool, and "site" when you change the home page.
   Keys match the folder names. Everything is at 0.2.0 as of 25 Nov 2026. */
window.CP_VERSIONS = {
  site: '0.2.0',
  'draw': '0.2.0',
  'dock': '0.2.0'
};
/* Display text: a number such as "0.2.0" is shown as "v0.2.0". */
window.cpVer = function (v) { return /^\d/.test(v) ? 'v' + v : v; };
