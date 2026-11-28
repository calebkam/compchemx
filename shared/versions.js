/* One version number for the whole platform. Any new feature or update bumps it (X.Y.Z, see the changelog). Individual tools carry no number of their own. */
window.CP_VERSIONS = { site: '0.5.1' };
/* Display text: a number such as "0.2.1" is shown as "v0.2.1". */
window.cpVer = function (v) { return /^\d/.test(v) ? 'v' + v : v; };
