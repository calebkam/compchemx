const D=require('../density.js');let fails=0;const ok=(c,m)=>{console.log((c?'ok   ':'FAIL ')+m);if(!c)fails++};
const rnd=(()=>{let s=7;return()=>{s=(s*16807)%2147483647;return s/2147483647}})();
/* a small planar ligand, a displaced copy and a far copy */
const L=[[10,10,10],[11.4,10,10],[12.1,11.2,10],[11.4,12.4,10],[10,12.4,10],[9.3,11.2,10],[13.5,11.2,10.3],[8,11.2,9.8]].map((p,i)=>({xyz:p,el:i%3?'C':'N'}));
const move=(a,d)=>a.map(x=>({el:x.el,xyz:[x.xyz[0]+d[0],x.xyz[1]+d[1],x.xyz[2]+d[2]]}));
const buf=D.synthMap(L,[40,40,40,90,90,90],0.3,0.7,0.4,rnd), g=D.parseMap(buf);
ok(g.m[0]===133&&g.dims[0]===133,'map header read ('+g.dims.join('x')+')');
const r0=D.rsr(g,L,L,{res:2,B:20}),r1=D.rsr(g,move(L,[0.4,0.3,0]),L,{res:2,B:20}),r2=D.rsr(g,move(L,[1.5,1.0,0]),L,{res:2,B:20}),r3=D.rsr(g,move(L,[8,0,0]),L,{res:2,B:20});
console.log('RSR exact',r0.rsr.toFixed(3),'0.5 A',r1.rsr.toFixed(3),'1.8 A',r2.rsr.toFixed(3),'far',r3.rsr.toFixed(3),'points',r0.n);
ok(r0.rsr<0.2,'exact pose fits the density well');
ok(r0.rsr<r1.rsr&&r1.rsr<r2.rsr&&r2.rsr<r3.rsr,'RSR rises as the pose moves away');
ok(r3.rsr>0.6,'a pose in empty space fits badly');
/* triclinic cell round trip */
const M=D.cellMatrix([30,35,40,80,95,100]),p=[3.2,-4.1,7.7];const f=D.toFrac(M,p),b=D.toCart(M,f);ok(Math.hypot(b[0]-p[0],b[1]-p[1],b[2]-p[2])<1e-9,'fractional/Cartesian round trip');
try{D.parseMap(new ArrayBuffer(2000));ok(false,'rejects garbage')}catch(e){ok(true,'rejects garbage: '+e.message)}
if(fails){console.log(fails+' failed');process.exit(1)}console.log('ALL PASSED');
