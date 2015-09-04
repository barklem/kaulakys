
pro test_v

v = findgen(100)*0.00001

crossmax = 0.095 * 2900 * !pi * 5.29177e-9^2 / (0.13^1.5) / sqrt(v) * 1e11
nmax = 1.31*sqrt(0.13/2./v)

plot, v, nmax

stop
end

pro test_v2

; if we combine equations for crossmax and nmax and eliminate velocity we get:

; this would imply nmax = 28 -> crossmax = 4, not 3 as seen in plot.
nmax = findgen(50)

crossmax = 0.103 * 2900 * !pi * 5.29177e-9^2 * 1e11 * nmax / 0.13^2.

plot, nmax, crossmax

stop
end




pro test_h
; compare with hoang-binh & vR 1995

T = 5000.
kT = T *1.38066e-23/4.35981e-18 ;  kT in au  = T in au (k=1)
A = 24.3  ; 
mt = A
mp = 1.008
mu = mp*mt/(mt+mp) * 1822.89
v = sqrt(8.*kT/!pi/mu) 
print, 'v = ', v  
E = 0.5*mu*v*v   ; corresponding CofM energy



L_singlet = 5.965     ; scattering lengths in au, from Schwartz 1961 (agree with Bhatia 2007)
L_triplet = 1.7686
L_eff = sqrt(0.25 * (L_singlet^2.) + 0.75 * (L_triplet^2.))

; this is 3.3528 and agrees with quoted value of 3.35 in Hoang-Binh


print, '7i -> 7l'

Elim = 61671.
E1 = 59430.517  ; 7i
E2 = 59428.853  ; 7h
dE = (E1-E2) / 8065. / 27.2116; in au


nstarupp = sqrt(109678./(Elim-E1))
nlow = 7.
llow = 5.
lupp = 6.

;goto, six

def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.

y = dE*nlow/v

;;;; need initial state -> final

inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = inel4 / 1510.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' ; broadening cross section

llow = 4.
E2 = 59423.537 ; 7g
dE = (E1-E2) / 8065. / 27.2; in au
def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.
y = dE*nlow/v
inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = inel4 / 1130.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' ; broadening cross section

llow = 3.
E2 = 59400. ; 7f
dE = (E1-E2) / 8065. / 27.2; in au
def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.
y = dE*nlow/v
inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = inel4 / 847.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' ; broadening cross section

llow = 2.
E2 = 59318. ; 7d
dE = (E1-E2) / 8065. / 27.2; in au
nstarlow = sqrt(109678./(Elim-E2))
def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.
y = dE*nlow/v
inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = inel4 / 546.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' ; broadening cross section

llow = 1.
E2 = 58476. ; 7p
dE = (E1-E2) / 8065. / 27.2; in au
def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.
y = dE*nlow/v
inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = 0.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' ; broadening cross section

llow = 0.
E2 = 57855.214 ; 7s
dE = (E1-E2) / 8065. / 27.2; in au
def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.
y = dE*nlow/v
inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = 0.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' ; broadening cross section

print, 'Omont      = ', 2.*!pi/7.^3.*L_eff^2. / v^2.

six:

print, '6h -> 6l'

Elim = 61671.
E1 = 58618.942  ; 6h
E2 = 58610.795  ; 6g
dE = (E1-E2) / 8065. / 27.2; in au

nstarupp = sqrt(109678./(Elim-E2))
nlow = 6.
llow = 4.
lupp = 5.

def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.

y = dE*nlow/v

inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = inel4 / 2410.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' 


llow = 3.
E2 = 58575.527 ; 6f
dE = (E1-E2) / 8065. / 27.2; in au
def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.
y = dE*nlow/v
inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = inel4 / 1770.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' 

llow = 2.
E2 = 58442.835 ; 6d
dE = (E1-E2) / 8065. / 27.2; in au
def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.
y = dE*nlow/v
inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = inel4 / 865.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' 

llow = 1.
E2 = 57017.078 ; 6p
dE = (E1-E2) / 8065. / 27.2; in au
def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.
y = dE*nlow/v
inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = 0.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' 

llow = 0.
E2 = 55891.80 ; 6s
dE = (E1-E2) / 8065. / 27.2; in au
def1 = 33.6/4./nlow/nlow*(3*nlow*nlow-llow*(llow+1))/((llow+1.5)*(llow+1.)*(llow+0.5)*llow*(llow-0.5))
def2 = 33.6/4./nlow/nlow*(3*nlow*nlow-lupp*(lupp+1))/((lupp+1.5)*(lupp+1.)*(lupp+0.5)*lupp*(lupp-0.5))
dE2 = abs(def1-def2)/nlow^3.
y = dE*nlow/v
inel = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1)
inel2 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, 0., E, cofm=1)
inel3 = kaulakys_crossh(A, nstarupp, lupp, nlow, llow, dE, E, cofm=1, version=2)
inel4 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE, E, cofm=1)
inel5 = kaulakys_crossh_best(A, nstarupp, nlow, lupp, nlow, llow, dE2, E, cofm=1)
ratio1 = 0.
print, llow, dE, dE2, y, inel, inel2, inel3, inel4, inel5, ratio1, format = '(i3, 3e13.5, 5f10.0, f10.3)' 

print, 'Omont      = ', 2.*!pi/6.^3.*L_eff^2. / v^2.

stop
end

pro test_kaulr


t2 = systime(1)
T = [2000., 4000., 6000., 8000., 10000]
dE = 1. ; eV
llow = 0
lupp = 1
nupp = 11
nstarlow = 7.9
mp = 1.
mt = 40.

E = T *1.38e-23/4.35981e-18 ;  energy in au  = T in au
mu = mp*mt/(mt+mp) * 1836.15
v = sqrt(2.*E/mu)
L_singlet = 5.965     ; scattering lengths in au, from Schwartz 1961 (agree with Bhatia 2007)
L_triplet = 1.7686
L_eff = sqrt(0.25 * (L_singlet^2.) + 0.75 * (L_triplet^2.))

r1 = kaulakys_rateh(mt, nstarlow, llow, nupp, lupp, dE, T)
r2 = kaulakys_cross86(mt, nstarlow, 0, nupp, mp, L_eff, dE/27.2, E, cofm=1) * v
r3 = kaulakys_cross86b(mt, nstarlow, 0, nupp, mp, L_eff, dE/27.2, E, cofm=1) * v 

print, r1
print, r2 * 5.29e-9^3/2.418e-17
print, r3* 5.29e-9^3/2.418e-17 ; in cgs
; not totally good comparison since r2 and r3 include all l values of nupp
print, systime(1) - t2, ' SECONDS'
stop
end

pro test_kaul2

t2 = systime(1)

T = 400000.  ; in K
mt = 86.
mp = 86.
lscat = 26.9

E = T *1.38e-23/4.35981e-18 ;  energy in au  = T in au
mu = mp*mt/(mt+mp) * 1836.15
v = sqrt(2.*E/mu)
nstar = 5.1
dE = 100.13/(5^3.)
cross = kaulakys_cross(85, nstar, 0, fix(nstar), 1, 85, 26.9, dE, E, cofm=1)
print, nstar, 1, cross* (5.29177e-9^2.) *1e11
print, systime(1) - t2, ' SECONDS'
end

pro test_kaul
; try to reproduce figure 2 of 1986 paper

T = 400.  ; in K
mt = 86.
mp = 86.
lscat = 26.9   ; L = sqrt(2900.*pi/4/pi)

; we need the collision energy corresponding to the average relative velocity

kT = T *1.38d-23/4.35981e-18 ;  energy in au  = T in au, which is energy of most probable velocity
mu = mp*mt/(mt+mp) * 1.66e-27 / 9.11e-31 ; reduced mass in au
vprob = sqrt(2.*kT/mu)
vbar = sqrt(8/!pi*kT/mu)
v = vbar
;v = vprob
E = 0.5*mu*v*v   

; I get good agreement with predicted maximum position and magnitude.
; smallish difference in plot though...
; note we can't solve this trivially
; if we adjust v to achieve nmax=28 as in figure, then cross max gets even larger too,
; and seems smaller in figure.  Delta n* is only other parameter and is 0.13

crossmax = 0.095 * 2900 * !pi * 5.29177e-9^2 / (0.13^1.5) / sqrt(v) * 1e11
nmax = 1.31*sqrt(0.13/2./v)

print, ' E = ', E
print, ' v = ', v
print, ' nstarmax = ', nmax
print, ' crossmax = ', crossmax

;stop

nlow = indgen(60) + 3 - 0.13
n = n_elements(nlow)
cross = nlow*0.d0
cross2 = nlow*0.d0
cross3 = nlow*0.d0

set_plot, 'ps'
device, file = 'idl.eps'
for i = 0, n-1 do begin
   nupp = round(nlow[i])*1.
    
   dE = 0.13/(nlow[i]^3.)
   for lupp = 1, nupp-1 do begin  ; other l will have different defects
       cross[i] = cross[i] + kaulakys_cross(mt, nlow[i], 0, nupp, lupp, mp, lscat, dE, E, cofm=1)
       ;print, nlow[i], lupp, cross[i]* (5.29177e-9^2.) *1e11
   endfor
   cross2[i] = kaulakys_cross86(mt, nlow[i], 0, nupp, mp, lscat, dE, E, cofm=1)
   cross3[i] = kaulakys_cross86b(mt, nlow[i], 0, nupp, mp, lscat, dE, E, cofm=1)
   print, nlow[i], nupp, cross[i]* (5.29177e-9^2.) *1e11, cross2[i]* (5.29177e-9^2.) *1e11, cross3[i]* (5.29177e-9^2.) *1e11
endfor


plot, nlow, cross * (5.29177e-9^2.) *1e11, charsize = 2
oplot, nlow, cross2 * (5.29177e-9^2.) *1e11, linestyle = 2
oplot, nlow, cross3 * (5.29177e-9^2.) *1e11, linestyle = 5, thick = 2
oplot, [nmax], [crossmax], psym = 4, symsize = 3

device, /close_file
set_plot, 'x'

stop
nn = 100
eta = dindgen(1000)*0.01+0.01
pt = eta / nn
y = pt * 0.d0
for i = 0, 999 do y[i] = Inl(nn,0,pt[i]) 
;plot, eta, y, /ylog, /xlog
;oplot, eta, 2./!pi * (atan(1./eta) - eta * alog(1.+1./(eta*eta))), linestyle = 2


nn = 5
pp = dindgen(100)*0.01+1./(2.*nn*nn)
;plot, pp, 4*nn/!pi/(pp*pp*(1.+nn*nn*pp*pp)^2.), linestyle = 2, charsize = 2
ff = Fnl(mt,0,nn,0,pp)
;oplot, pp, ff*ff

stop
end

function hyper2F1old, a, b, c, z
; calculates the 2F1 hypergeometric series
; this is an older version, non-vectorised version, but well tested.

if (a eq 0.) or (b eq 0.) or (c eq 0.) then return, 0.d0

errlim = 1d-16  
err = 1d10
sum = 0.d0
n = 0l
while (err gt errlim) do begin
  prod = 1.d0 
  ; term by term evaluation of product for numerical precision  
  if n ge 1 then begin
     iarr = dindgen(n)
     prod =  product(1.d0 / (iarr+1.) * (a+iarr) / (c+iarr) * (b+iarr))
  endif
  term = prod * double(z)^n
  if n ge 1 then err = abs(term/sum)
  sum = sum + term
  n = n + 1
endwhile

return, sum
end


pro testhankel
q = dindgen(100)/100.+1
nu = 10.0
l = 9
t = 4
!p.multi=[0,1,3]
y1 = q*0.d0
for i = 0, n_elements(q)-1 do y1[i] = hankel(q[i],nu,l,t)
y2 = sqrt(2./!pi) * gamma(nu-t+1) * ((1./nu/nu + q*q)^((-nu+t-1)/2.))*sin((nu-t+1)*atan(q*nu))
plot, q, y1, charsize = 2
plot, q, y2, psym = 5, charsize = 2
plot, q, y1/y2, charsize = 2, thicK = 3
!p.multi=0
stop
end

pro testhyper

z = dindgen(100)/100.
a = .5
b = 1.
c = 1.5
h1 = z*0.d0
for i = 0, n_elements(z)-1 do h1[i] = z[i] * hyper2F1(a,b,c,z[i]*z[i])
h2 = 0.5*alog((1+z)/(1-z))
plot, z, h2, linestyle=2, thick = 3
oplot, z, h1
stop
end

pro testgwf
t = systime(1)
q = (dindgen(100))/40
window, 0, xsize = 1900
!p.multi = [0,5,1]
!p.charsize = 3
; Na 4s
ns = 2.64
plot, q, q*gwf(ns, 0, q), ytitle = "q*g(q)", xtitle = 'q (au)', title = 'Na (4s)'
oplot, q, q*gwf(2.0, 0, q), linestyle = 2
oplot, q, q*gwf(3.0, 0, q), linestyle = 2
;oplot, q, (q*gwf(2.2, 0, q))^2, linestyle = 3
;oplot, q, (q*gwf(2.8, 0, q))^2, linestyle = 3
oplot, q, q*Fnl(1., 0, 3, 0, q), psym = 4
; Na 4p
ns = 3.13
plot, q, q*gwf(ns, 1, q), xtitle = 'q (au)', title = 'Na (4p)'
oplot, q, q*gwf(3.0, 1, q), linestyle = 2
oplot, q, -q*Fnl(1., 0, 3, 1, q), psym = 4
; Mg 6s
ns = 4.358
plot, q, q*gwf(ns, 0, q), title = 'Mg (6s)'
oplot, q, q*gwf(4.0, 2, q), linestyle = 2
oplot, q, -q*Fnl(1., 0., 4.0, 2, q), psym = 4
; Mg 6p
ns = 4.86
plot, q, q*gwf(ns, 1, q), title = 'Mg (6p)';, xr=[0,1.0], /xs
oplot, q, q*gwf(5.0, 1, q), linestyle = 2
oplot, q, -q*Fnl(1., 0, 5, 1, q), psym = 4
; Cs 10p
ns = 6.43
plot, q, q*gwf(ns, 1, q), title = 'Cs (10p)';, xr=[0,1.0], /xs
oplot, q, q*gwf(6.0, 1, q), linestyle = 2
oplot, q, q*Fnl(1., 0, 6, 1, q), psym = 4
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwfl
t = systime(1)
q = [(dindgen(20)+1), 50.]
;window, 0, xsize = 1700, ysize = 1000
set_plot, 'ps'
device, file = 'wfs.eps', ysize = 30, xsize = 45
!p.multi = [0,6,4]
!p.charsize = 2
; Na 4s
l = 0
ns = 2.1
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 5.7
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs 
print, systime(1) - t, ' SECONDS'
ns = 12.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'

l = 2
ns = 5.7
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 10.0
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS' 
ns = 12.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'


l = 4
ns = 5.7
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 10.0
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS' 
ns = 12.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'

l = 7
ns = 8.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 10.0
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS' 
ns = 12.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'

l = 13
ns = 14.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
ns = 15.0
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs 
ns = 16.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
;plot, q, (q*gwf(2.2, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
;oplot, q, (q*gwf(2.8, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
;oplot, q, (q*gwf(2.0, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
;oplot, q, (q*gwf(3.0, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
;oplot, q, (q*Fnl(1., 0, 3, 0, q))^2, psym = 4
print, systime(1) - t, ' SECONDS'
device, /close_file
set_plot, 'x'
stop
; Na 4p
ns = 3.13
plot, q, (gwf(ns, 1, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
oplot, q, (gwf(3.0, 1, q))^2, linestyle = 2
oplot, q, (Fnl(1., 0, 3, 1, q))^2, psym = 4
; Mg 6s
ns = 4.358
plot, q, (gwf(ns, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
oplot, q, (gwf(4.0, 0, q))^2, linestyle = 2
oplot, q, (Fnl(1., 0, 4, 0, q))^2, psym = 4
; Mg 6p
ns = 4.86
plot, q, (gwf(ns, 1, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs ;, xr=[0,1.0], /xs
oplot, q, (gwf(5.0, 1, q))^2, linestyle = 2
oplot, q, (Fnl(1., 0, 5, 1, q))^2, psym = 4
; Cs 10p
ns = 6.43
plot, q, (gwf(ns, 1, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs ;, xr=[0,1.0], /xs
oplot, q, (gwf(6.0, 1, q))^2, linestyle = 2
oplot, q, (Fnl(1., 0, 6, 1, q))^2, psym = 4
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwfl2
t = systime(1)
pmin = 1e-4
pmax = 1000.
npts = 10000
pstep = (alog(pmax) - alog(pmin))/(npts-1)
q = exp(alog(pmin) + pstep*dindgen(npts))
set_plot, 'ps'
device, file = 'wfs.eps', xsize = 28, xoffset=1, yoffset = 29, /landscape
!p.multi = 0
!p.charsize = 1.5

ns = 20.0
yy = gwf(ns, 0, q)
plot, q, (yy)^2, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')
for ns = 3.0, 20.0, 1.0 do begin
yy = gwf(ns, 0, q)
oplot, q, (yy)^2, linestyle = 5
endfor

for ns = 3.0, 20.0, 1.0 do begin

yy = gwf(ns, 0, q)
plot, q, (yy)^2, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')

for l = 1, min([13, fix(ns)-1]), 1 do begin
oplot, q, gwf(ns, l, q)^2, linestyle = 5
end

for l = 0, min([13, fix(ns)-1]), 1 do begin

yy = gwf(ns, l, q)
plot, q, (yy)^2, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
oplot, q, gwf(ns-0.5, l, q)^2, linestyle = 5
oplot, q, gwf(ns-0.4, l, q)^2,linestyle = 5
oplot, q, gwf(ns-0.3, l, q)^2, linestyle = 5
oplot, q, gwf(ns-0.2, l, q)^2, linestyle = 5
oplot, q, gwf(ns-0.1, l, q)^2, linestyle = 5
oplot, q, gwf(ns+0.5, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.4, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.3, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.2, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.1, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.01, l, q)^2, linestyle = 1, thick = 10

oplot, [2/!Pi/ns,2/!Pi/ns], [1e-30,1e30] 
oplot, [2/!Pi/ns,2/!Pi/ns]*2, [1e-30,1e30], linestyle = 5 
oplot, [2/!Pi/ns,2/!Pi/ns]*4, [1e-30,1e30], linestyle = 1 
oplot, [2/!Pi/ns,2/!Pi/ns]*.5, [1e-30,1e30], linestyle = 5 
oplot, [2/!Pi/ns,2/!Pi/ns]*.25, [1e-30,1e30], linestyle = 1 
print, systime(1) - t, ' SECONDS'
end
end

device, /close_file
set_plot, 'x'
!p.multi = 0
!p.charsize = 1
end

pro testgwfl3
t = systime(1)
pmin = 1e-4
pmax = 1000.
npts = 10000
pstep = (alog(pmax) - alog(pmin))/(npts-1)
q = exp(alog(pmin) + pstep*dindgen(npts))
set_plot, 'ps'
device, file = 'wfs_new.eps', xsize = 28, xoffset=1, yoffset = 29, /landscape
!p.multi = 0
!p.charsize = 1.5

ns = 20.0
yy = gwfsq(ns, 0, q)
plot, q, yy, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')
for ns = 3.0, 20.0, 1.0 do begin
yy = gwfsq(ns, 0, q)
oplot, q, yy, linestyle = 5
endfor

for ns = 3.0, 20.0, 1.0 do begin

yy = gwfsq(ns, 0, q)
plot, q, yy, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')

for l = 0, min([13, fix(ns)-1]), 1 do begin
oplot, q, gwfsq(ns, l, q), linestyle = 5
end

for l = 0, min([13, fix(ns)-1]), 1 do begin

print, ns, l

yy = gwfsq(ns, l, q)
plot, q, yy, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
oplot, q, gwfsq(ns-0.5, l, q), linestyle = 5
oplot, q, gwfsq(ns-0.4, l, q),linestyle = 5
oplot, q, gwfsq(ns-0.3, l, q), linestyle = 5
oplot, q, gwfsq(ns-0.2, l, q), linestyle = 5
oplot, q, gwfsq(ns-0.1, l, q), linestyle = 5
oplot, q, gwfsq(ns+0.5, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.4, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.3, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.2, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.1, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.01, l, q), linestyle = 1, thick = 10

oplot, [2/!Pi/ns,2/!Pi/ns], [1e-30,1e30] 
oplot, [2/!Pi/ns,2/!Pi/ns]*2, [1e-30,1e30], linestyle = 5 
oplot, [2/!Pi/ns,2/!Pi/ns]*4, [1e-30,1e30], linestyle = 1 
oplot, [2/!Pi/ns,2/!Pi/ns]*.5, [1e-30,1e30], linestyle = 5 
oplot, [2/!Pi/ns,2/!Pi/ns]*.25, [1e-30,1e30], linestyle = 1 
print, systime(1) - t, ' SECONDS'
end
end

device, /close_file
set_plot, 'x'
!p.multi = 0
!p.charsize = 1
end


pro testgwf2
t = systime(1)
q = (dindgen(1000))/10000.
window, 0, xsize = 1900
!p.multi = [0,5,1]
!p.charsize = 3
ns = 20.
plot, q, q*gwf(ns, 0, q), ytitle = "q*g(q)", xtitle = 'q (au)', title = '20s'
ns = 20.
plot, q, q*gwf(ns, 1, q), xtitle = 'q (au)', title = '20p'
ns = 25.
plot, q, q*gwf(ns, 0, q), xtitle = 'q (au)', title = '25s'
ns = 25.
plot, q, q*gwf(ns, 1, q), xtitle = 'q (au)', title = '25p';, xr=[0,1.0], /xs
ns = 30.
plot, q, q*gwf(ns, 1, q), xtitle = 'q (au)', title = '30p';, xr=[0,1.0], /xs
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwf3
t = systime(1)
print, systime()
q = (dindgen(10))/10000. + 0.0495001
q2 = (dindgen(1000))/1000000. + 0.0495
window, 0, xsize = 1900
!p.multi = [0,2,1]
!p.charsize = 3
ns = 10.
y1 = q*gwf(ns, 0, q)
y1a = -q2*Fnl(1., 0, 10, 0, q2)
plot, q, y1, ytitle = "q*g(q)", xtitle = 'q (au)', title = '20s'
oplot, q, y1a, psym = 4
print, systime(1) - t, ' SECONDS'
y2 = q*gwf(ns, 1, q)
y2a = q2*Fnl(1., 0, 20, 1, q2)
plot, q, y2, xtitle = 'q (au)', title = '20p'
oplot, q2, y2a, psym = 4
print, systime(1) - t, ' SECONDS'
print, systime()
!p.multi = 0
!p.charsize = 1
end


pro testgwf4
t = systime(1)
q = (dindgen(200))/1000.
window, 0, xsize = 1900
!p.multi = [0,5,1]
!p.charsize = 3
plot, q, q*gwf(15, 0, q), ytitle = "q*g(q)", xtitle = 'q (au)'
oplot, q, q*Fnl(1., 0, 15, 0, q), psym = 4

plot, q, q*gwf(20, 0, q)
oplot, q, -q*Fnl(1., 0, 20, 0, q), psym = 4

plot, q, q*gwf(25, 0, q)
oplot, q, q*Fnl(1., 0, 25, 0, q), psym = 4

plot, q, q*gwf(30, 0, q)
oplot, q, -q*Fnl(1., 0, 30, 0, q), psym = 4

plot, q, q*gwf(50, 0, q)
oplot, q, q*Fnl(1., 0, 50, 0, q), psym = 4

print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwf5
t = systime(1)
q = (dindgen(1000))/2000.
window, 0, xsize = 1900
!p.multi = [0,1,1]
!p.charsize = 3
plot, q, q*gwf(9, 0, q), yr = [-4,4]
oplot, q, q*Fnl(1., 0, 9, 0, q), psym = 4
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwf6
t = systime(1)
q = (dindgen(500))/5000.
window, 0, xsize = 1900
!p.multi = [0,1,1]
!p.charsize = 3
plot, q, q*Fnl(1., 0, 30, 0, q)
;oplot, q, q*Fnl(1., 0, 30, 0, q), linestyle = 1
oplot, q, -q*gwf(30., 0, q), linestyle = 5
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end
