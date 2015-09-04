pro test

T = [2000, 3000, 4000, 5000, 6000, 7000, 8000]
;T = [2000, 4000, 6000, 8000]
ts = {A:24, N:3, L:2, ND:4, LD:1, NSTAR:2.6809385, DE:0.17954104}


print, t
tt = systime(1)
met1 = kaulakys_rateh (T, ts, method=1)
tt = systime(1)-tt
print, ' Method 1 : ', tt
print, met1
tt = systime(1) 
met2a = kaulakys_rateh (T, ts, method=2, npts=10, plt=1)
tt = systime(1)-tt
print, ' Method 2 (10) : ', tt
print, met2a
tt = systime(1) 
met2 = kaulakys_rateh (T, ts, method=2, npts=50, plt=1)
tt = systime(1)-tt
print, ' Method 2 (50) : ', tt
print, met2


plot, t, met1, linestyle = 1, /ylog
oplot, t, met2, linestyle = 2
oplot, t, met2a, linestyle = 2

tt = systime(1) 
met3 = kaulakys_rateh (T, ts, method=3)
tt = systime(1)-tt
print, ' Method 3 : ', tt
print, met3

oplot, t, met3, linestyle = 3

end