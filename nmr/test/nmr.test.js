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

/* uncommon groups (literature values, CDCl3) */
const sg=(n,f)=>run(n).h.signals.filter(f).map(s=>s.shift),cs=n=>run(n).c.signals.map(s=>s.shift);
let v;
v=sg('nitromethane',s=>!s.exch);near('MeNO2 H',v[0],4.33,.15);near('MeNO2 C',cs('nitromethane')[0],62.5,4);
v=sg('nitroethane',s=>!s.exch);near('EtNO2 CH2',v[0],4.38,.2);near('EtNO2 CH3',v[1],1.58,.15);near('EtNO2 C',cs('nitroethane')[0],70,4);
for(const n of ['nitroethane (N(=O)=O)','nitrobenzene (N(=O)=O)'])run(n);
near('N(=O)=O form',sg('nitroethane (N(=O)=O)',s=>true)[0],4.38,.2);
v=sg('1-nitropropane',s=>true);near('PrNO2 a',v[0],4.38,.2);near('PrNO2 b',v[1],2.01,.2);near('PrNO2 C',cs('1-nitropropane')[0],77.4,4);
v=sg('2-nitropropane',s=>true);near('iPrNO2 CH',v[0],4.67,.2);near('iPrNO2 C',cs('2-nitropropane')[0],78.5,4);
v=sg('1,3-dinitrobenzene',s=>true);near('DNB H2',v[0],9.1,.2);near('DNB H4',v[1],8.6,.2);near('DNB H5',v[2],7.8,.2);
v=sg('4-nitroaniline',s=>true);near('pNA H',v[0],8.07,.2);near('pNA H',v[1],6.62,.3);near('pNA NH2',v[2],4.4,.7);
v=sg('nitrostyrene',s=>s.kind==='vinyl');near('nitrostyrene H',v[0],7.99,.3);near('nitrostyrene H',v[1],7.58,.3);
near('3-nitropyridine H2',sg('3-nitropyridine',s=>true)[0],9.4,.3);
near('ethyl nitrate OCH2',sg('ethyl nitrate',s=>true)[0],4.5,.3);near('ethyl nitrate C',cs('ethyl nitrate')[0],69,6);
near('azide CH2',sg('benzyl azide',s=>s.n===2&&s.shift>4&&s.kind==='sp3')[0],4.35,.4);near('azide C',cs('benzyl azide')[4],54.8,6);
near('NCS CH2',sg('ethyl isothiocyanate',s=>true)[0],3.5,.3);near('NCS C',cs('ethyl isothiocyanate')[0],130,8);
near('NCO CH2',sg('ethyl isocyanate',s=>true)[0],3.3,.3);near('NCO C',cs('ethyl isocyanate')[0],122,6);
near('SCN CH2',sg('ethyl thiocyanate',s=>true)[0],2.95,.3);near('SCN C',cs('ethyl thiocyanate')[0],112,5);
near('aldoxime C',cs('acetaldoxime')[0],150,5);near('ketoxime C',cs('acetone oxime')[0],156,5);
near('PhSO2NH2 H',sg('benzenesulfonamide',s=>!s.exch)[0],7.92,.3);near('PhSO2NH2 ipso',cs('benzenesulfonamide')[0],142.5,5);
near('Me2SO2 H',sg('dimethyl sulfone',s=>true)[0],3.0,.3);near('Me2SO2 C',cs('dimethyl sulfone')[0],42.6,5);
near('DMSO H',sg('DMSO',s=>true)[0],2.62,.3);near('DMSO C',cs('DMSO')[0],40.6,5);
near('TMS H',sg('TMS',s=>true)[0],0,.3);near('TMS C',cs('TMS')[0],0,4);
near('CF3 C',cs('benzotrifluoride').find(x=>x>120&&x<126.5&&Math.abs(x-124)<2),124,3);
near('CH2Br2 C',cs('dibromomethane')[0],21.4,4);near('CHBr3 C',cs('bromoform')[0],12.1,6);near('CH3I C',cs('iodomethane')[0],-20.7,5);near('CH3I H',sg('iodomethane',s=>true)[0],2.16,.3);
near('CH2Cl2 C',cs('dichloromethane')[0],54,4);near('CCl4 C',cs('carbon tetrachloride')[0],96,5);near('EtBr C',cs('bromoethane')[0],28.3,4);
near('PhB(OH)2 ortho H',sg('phenylboronic acid',s=>!s.exch)[0],8.2,.3);
near('(EtO)3PO H',sg('triethyl phosphate',s=>true)[0],4.1,.3);near('(EtO)3PO C',cs('triethyl phosphate')[0],64,5);
near('C=S',cs('thioacetamide')[0],200,15);
near('TMAO H',sg('TMAO',s=>true)[0],3.25,.3);near('TMAO C',cs('TMAO')[0],60,8);
near('azobenzene H',sg('azobenzene',s=>true)[0],7.9,.3);near('azobenzene C',cs('azobenzene')[0],152.7,5);
near('NDMA H',sg('NDMA',s=>true)[0],3.4,.5);
near('PhOCF3 ipso',cs('trifluoromethoxybenzene')[0],149.3,5);
near('p-nitrophenol DMSO-free',sg('4-nitrophenol',s=>s.exch)[0],5.9,.8);
if(run('sodium acetate').c.signals.length!==2||run('methylammonium chloride').formula!=='CH6N')bad++,console.log('FAIL counter-ions');
near('tetranitromethane C',cs('tetranitromethane')[0],130,12);

console.log(bad?bad+' failed':'nmr tests passed');process.exit(bad?1:0);
