const N=require('../nmr.js'),M=require('./molfiles.json');let bad=0;
const near=(n,got,exp,tol)=>{const ok=Math.abs(got-exp)<=tol;if(!ok){bad++;console.log('FAIL',n,got.toFixed(2),'expected',exp,'+/-',tol)}return ok};
function run(name){const r=N.predict(M[name]);return r}
const show=process.argv[2]==='-v';
for(const k of Object.keys(M)){const r=run(k);if(show){console.log('\n'+k,r.formula);r.h.signals.forEach(s=>console.log('  H',s.shift.toFixed(2),s.n+'H',s.mult,s.label));console.log('  C',r.c.signals.map(s=>s.shift.toFixed(0)+'('+s.n+')').join(' '))}}
const H=(n)=>run(n).h.signals.filter(s=>!s.exch).map(s=>[s.shift,s.n,s.mult]),C=n=>run(n).c.signals.map(s=>s.shift);
let h=H('ethanol');near('EtOH CH3',h[1][0],1.2,.25);near('EtOH CH2',h[0][0],3.7,.3);
h=H('ethyl acetate');near('EtOAc OCH2',h[0][0],4.1,.2);near('EtOAc Ac',h[1][0],2.05,.2);near('EtOAc CH3',h[2][0],1.26,.2);
if(h[0][2]!=='q'||h[2][2]!=='t')bad++,console.log('FAIL EtOAc mult',h);
h=H('toluene');near('tol CH3',h.find(x=>x[1]===3)[0],2.34,.2);
near('acetone',H('acetone')[0][0],2.17,.15);
h=H('benzaldehyde');near('PhCHO',h[0][0],10.0,.15);
near('CHCl3',H('chloroform')[0][0],7.26,.4);
{const v=run('styrene').h.signals.filter(x=>x.kind==='vinyl').map(x=>x.shift);near('styrene CH=',v[0],6.72,.3);near('styrene cis',v[1],5.75,.3);near('styrene trans',v[2],5.25,.3);}
h=H('pyridine');near('py H2',h[0][0],8.6,.2);near('py H4',h[1][0],7.64,.25);near('py H3',h[2][0],7.25,.2);
h=H('naphthalene');near('naph a',h[0][0],7.84,.15);near('naph b',h[1][0],7.48,.15);
h=H('furan');near('furan a',h[0][0],7.44,.15);near('furan b',h[1][0],6.38,.15);
h=H('p-xylene');if(h.length!==2||h[0][2]!=='s'&&h[0][1]!==4)console.log('p-xylene',h);
let c=C('toluene');near('tol C ipso',c[0],137.8,3);near('tol CH3',c[c.length-1],21.4,3);
c=C('acetone');near('acetone CO',c[0],206.7,4);near('acetone CH3',c[1],30.8,4);
c=C('cyclohexane');near('cyclohexane',c[0],27,2);
c=C('ethanol');near('EtOH CH2',c[0],58,4);near('EtOH CH3',c[1],18,4);
c=C('pyridine');near('py C2',c[0],149.9,2);
c=C('acetonitrile');near('CN',c[0],118,2);
c=C('benzoic acid');near('COOH',c[0],172.6,5);
console.log(bad?bad+' failed':'nmr tests passed');process.exit(bad?1:0);
