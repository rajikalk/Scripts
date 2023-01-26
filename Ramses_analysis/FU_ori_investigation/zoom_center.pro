restore, '/lustre/astro/troels/IMF_512/stars_red_512.sav'
is=196
ii = where(s[is,*].m gt 0) & ss=reform(s[is,ii[0]])
print, (ss.t - ss.time)*[ss.ux,ss.uy,ss.uz]+[ss.x,ss.y,ss.z]
print, (ss.t - ss.time)
tsnap=0.110059438406200E+01
tmid = 0.5*(tsnap + ss.t)
dt = tmid - ss.time
print, dt*[ss.ux,ss.uy,ss.uz]+[ss.x,ss.y,ss.z]
print, tsnap, tmid, ss.t, ss.time, dt, ss.t - ss.time

END
