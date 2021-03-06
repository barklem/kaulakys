; routines for calculation of Kaulakys free electron model
;
; Paul Barklem
; written Jan-March 2010
; lots of new parts Jan 2013
;
; Notes: 1. various parameters in this code may need tweaking for your particular case (esp. pmax and npts in Inl)
;

function scatamp_h, p, x
; gives the SQUARE of the scattering amplitude f_e(p,x) in au for e+H collisions 
; a function of momentum p in au, and of x= cos theta
;
;
;  Notes, forcing A > 0 doesn't change anything in tests (will go < 0 only at very large p)
;  Spline goes crazy;) 

ptab = [0., 0.1, 0.2, 0.3]
Atab_singlet = [35.64, 30.58, 19.36, 10.96]
Btab_singlet = [0., -1.93, -0.88, -0.14]      ; s-p coupling
Atab_triplet = [ 3.13,  4.07,  4.24,  3.98]
Btab_triplet = [0., -1.42, -2.66, -3.44]      ; s-p coupling

As = interpol(Atab_singlet, ptab, p)       ; linear interpolation
Bs = interpol(Btab_singlet, ptab, p)   

At = interpol(Atab_triplet, ptab, p)       ; linear interpolation
Bt = interpol(Btab_triplet, ptab, p)   

; forcing > 0 has no effect in tests, but is retained for safety
fsq_s = As + Bs * x > 0.
fsq_t = At + Bt * x > 0.

fsq = 0.25 * fsq_s + 0.75 * fsq_t

;stop
return, fsq
end

function ikern1, x, p=p, x1=x1, x2=x2
ans = x*0.
ind = where((1-x)*(x2-x)*(x-x1) gt 0., nind)
if nind gt 0 then ans[ind] = scatamp_h(p,x[ind]) / sqrt((1-x[ind])*(x2-x[ind])*(x[ind]-x1))
ind2 = where(~finite(ans), nind2)
if nind2 gt 0 then stop
return, ans
end

function ikern2, x, p=p, x1=x1, x2=x2
ans = x*0.
ind = where((1-x)*(-x1-x)*(x2+x) gt 0., nind)
if nind gt 0 then ans[ind] = scatamp_h(p,x[ind]) / sqrt((1-x[ind])*(-x1-x[ind])*(x2+x[ind]))
ind2 = where(~finite(ans), nind2)
if nind2 gt 0 then stop
return, ans
end


function igr1, p, x1, x2, xt2
ans = qpint1d('ikern1', x1, xt2, breakpoints = [1., x1, x2], functargs = {x1:x1, x2:x2, p:p})
return, ans
end

function igr2, p, x1, x2, xt1
ans = qpint1d('ikern2', -x2, xt1, breakpoints = [1., -x1, -x2], functargs = {x1:x1, x2:x2, p:p})
return, ans
end

function ksigma, p, pt, x1, x2
; p could be an array
n = n_elements(p)

xt1 = p*0. & xt2 = p*0.
i1  = p*0. & i2  = p*0.
sig = p*0. 

for i = 0, n-1 do begin
   xt1[i]= min([1-2*pt*pt/p[i]/p[i], -x1])
   xt2[i]= min([1-2*pt*pt/p[i]/p[i], x2])
   i1[i] = igr1(p[i], x1, x2, xt2[i])
   i2[i] = igr2(p[i], x1, x2, xt1[i])
   sig[i] = abs(i1[i]) + abs(i2[i])
   ;if (i1[i] lt 0) or (i2[i] lt 0) then stop
   if (~finite(i1[i])) or (~finite(i2[i])) then stop
endfor

return, sig
end

function kkernel, p, pt=pt, x1=x1, x2=x2, nstar=nstar, l=l
return, ksigma(p, pt, x1, x2) * gwfsq(nstar,l,p) * p*p
end

function kaulakys_crossh_best, A, nstar, n, l, nd, ld, dE, E, cofm=cofm, version=version

APert = 1.008
mu = A*APert*1836.15/(A+APert)                   ; reduced mass for target-pertuber
if keyword_set(cofm) then begin
   Ecm = E
endif else begin
   Ecm = double(A)/(A+APert) * E
endelse

; allow array of energies
cross = dblarr(n_elements(Ecm))
for i = 0, n_elements(Ecm)-1 do begin
  v = sqrt(2.d0*Ecm[i]/mu)                         ; relative velocity in au
  pt = dE / (2.d0*v)                               ; momentum transfer in au 
  if Ecm[i] gt dE then begin
     
   x1 = (sqrt(n*n-(l+0.5d0)^2.)*sqrt(nd*nd-(ld+0.5)^2.)-(l+0.5)*(ld+0.5))/(n*nd)
   x2 = (sqrt(n*n-(l+0.5d0)^2.)*sqrt(nd*nd-(ld+0.5)^2.)+(l+0.5)*(ld+0.5))/(n*nd)
   
   inf = !values.d_infinity
   igral = qpint1d('kkernel', pt, +inf, functargs = {pt:pt, x1:x1, x2:x2, nstar:nstar, l:l})
   cross[i] = (2.*ld+1.)*sqrt(2.)/2./nd^5./v^2. * igral
   
   
   if (~finite(cross[i])) then begin
   npts= 1000 & pmin = pt & pmax = 2.
   p = dindgen(npts) * (pmax-pmin)/(npts-1) + pmin
   ii = kkernel(p, pt=pt, x1=x1, x2=x2, nstar=nstar, l=l)
   stop
   endif
   
  endif else begin
     cross[i] = 0.
  endelse
endfor

return, cross
end

function Inl_kernel, p
; calculates the kernel of the Inl integral 

common QUANTUM, nstarp, lp
common SHARE, pt

wfsq = gwfsq(nstarp,lp,p)
kernel = wfsq*(p-pt)*p

return, kernel
end

function Inl, nstar, l, pt
; calculates the Inl function, a type of overlap integral

common QUANTUM, nstarp, lp
nstarp = nstar
lp = l

pmsq = 2.d0/!pi/nstarp                    ; <p>
pmin = max([double(pt), 1e-6])
pmax = max([11., pt*25]) 

npts = 10000 ; 20000 gave no significant improvement for this pmax choice at high n
if pt lt 0.0001 then npts = 20000            
; log grid - makes little diff wrt linear grid, npts and pmax are most important
pstep = (alog(pmax) - alog(pmin))/(npts-1)
p = exp(alog(pmin) + pstep*dindgen(npts))
; linear grid
;p = dindgen(npts) * (pmax-pmin)/(npts-1) + pmin
kern = Inl_kernel(p)
igral = int_tabulated(p, abs(kern), /double)
;plot, p[1:npts-1], kern[1:npts-1], /ylog, /xlog, title = string(nstar, '(f6.1)') + string(pt, '(e16.4)')+ string(igral, '(e16.4)')
return, igral
end

function kaulakys_cross, A, nstar, l, nd, ld, APert, Lscat, dE, E, cofm=cofm 
; the cross section in scattering length approximation (eq. 18, 1991 paper)
; where:
; A is the atomic mass of target (e.g. 40 for Ca)
; n and l are quantum numbers for upper and lower states of the target
; nstar though, is the effective principal quantum number
; APert is the atomic mass of the perturber (e.g. 1 for H)
; Lscat is the scattering length of the perturber-electron system in au
; dE is the energy spacing of the levels in au
; E is the kinetic energy in the laboratory (target's) frame in au
;   if cofm is set then it's in the centre-of-mass frame

common SHARE, pt

mu = A*APert*1836.15/(A+APert)                   ; reduced mass for target-pertuber
if keyword_set(cofm) then begin
   Ecm = E
endif else begin
   Ecm = double(A)/(A+APert) * E
endelse

; allow array of energies
cross = dblarr(n_elements(Ecm))
for i = 0, n_elements(Ecm)-1 do begin
  v = sqrt(2.d0*Ecm[i]/mu)                         ; relative velocity in au
  pt = dE / (2.d0*v)                               ; momentum transfer in au 
  if Ecm[i] gt dE then begin
     cross[i] = 2.*!pi*Lscat^2.*(2.*ld+1.)/nd^5./v^2. * Inl(nstar, l, pt)
  endif else begin
     cross[i] = 0.
  endelse
endfor

return, cross
end

function kaulakys_cross2, A, nstar, l, nd, ld, APert, Lscat, dE, E, cofm=cofm 
; the cross section in scattering length approximation (eq. 12, 1991 paper)
; assume degeneracy
; where:
; A is the atomic mass of target (e.g. 40 for Ca)
; n and l are quantum numbers for upper and lower states of the target
; nstar though, is the effective principal quantum number
; APert is the atomic mass of the perturber (e.g. 1 for H)
; Lscat is the scattering length of the perturber-electron system in au
; dE is the energy spacing of the levels in au
; E is the kinetic energy in the laboratory (target's) frame in au
;   if cofm is set then it's in the centre-of-mass frame

common SHARE, pt

mu = A*APert*1836.15/(A+APert)                   ; reduced mass for target-pertuber
if keyword_set(cofm) then begin
   Ecm = E
endif else begin
   Ecm = double(A)/(A+APert) * E
endelse

; allow array of energies
cross = dblarr(n_elements(Ecm))
for i = 0, n_elements(Ecm)-1 do begin
  v = sqrt(2.d0*Ecm[i]/mu)                         ; relative velocity in au
  pt = 0. ;dE / (2.d0*v)                               ; momentum transfer in au 
  if Ecm[i] gt dE then begin
     ; n = nd here always
     x1 = (sqrt(nd*nd-(l+0.5)^2.)*sqrt(nd*nd-(ld+0.5)^2.)-(l+0.5)*(ld+0.5))/(nd*nd)
     x2 = (sqrt(nd*nd-(l+0.5)^2.)*sqrt(nd*nd-(ld+0.5)^2.)+(l+0.5)*(ld+0.5))/(nd*nd)
     k1 = sqrt((x2-x1)/(1-x1))
     k2 = sqrt((x2-x1)/(1+x2))
     Igral =  Inl(nstar, l, 0.)
     ell = ( ellipticK(k1)/sqrt(1-x1) + ellipticK(k2)/sqrt(1+x2) )
     
     cross[i] = Lscat^2.*(2.*ld+1.)*sqrt(2.)/nd^5./v^2. * igral * ell
     ;print, nstar, l, dE, v, pt, igral, ell   
     ;if ld eq 1 then stop   
  endif else begin
     cross[i] = 0.
  endelse
endfor

return, cross
end

function kaulakys_cross86, A, nstarlow, llow, nupp, APert, Lscat, dE, E, cofm=cofm 
; this is for nl -> n', from 1986 paper
; using the Inl calculated directly here

common SHARE, pt

mu = A*APert*1836.15/(A+APert)                   ; reduced mass for target-pertuber
if keyword_set(cofm) then begin
   Ecm = E
endif else begin
   Ecm = double(A)/(A+APert) * E
endelse

; allow array of energies
cross = dblarr(n_elements(Ecm))
for i = 0, n_elements(Ecm)-1 do begin
  v = sqrt(2.d0*Ecm[i]/mu)                          ; relative velocity in au
  pt = dE / (2.d0*v)                                ; momentum transfer in au 
  eta = pt * nstarlow
  if Ecm[i] gt dE then begin
     cross[i] = 2.*!pi*Lscat^2./nupp^3./v^2. * Inl(nstarlow, llow, pt)
  endif else begin
     cross[i] = 0.
  endelse
endfor

return, cross
end

function kaulakys_cross86b, A, nstarlow, llow, nupp, APert, Lscat, dE, E, cofm=cofm 
; this is for nl -> n', from 1986 paper
; using the analytic expression for Inl

common SHARE, pt

mu = A*APert*1836.15/(A+APert)                   ; reduced mass for target-pertuber
if keyword_set(cofm) then begin
   Ecm = E
endif else begin
   Ecm = double(A)/(A+APert) * E
endelse

; allow array of energies
cross = dblarr(n_elements(Ecm))
for i = 0, n_elements(Ecm)-1 do begin
  v = sqrt(2.d0*Ecm[i]/mu)                          ; relative velocity in au
  pt = dE / (2.d0*v)                                ; momentum transfer in au 
  eta = pt * nstarlow
  if Ecm[i] gt dE then begin
     cross[i] = 2.*!pi*Lscat^2./nupp^3./v^2. * 2./!pi * (atan(1./eta) - eta * alog(1.+1./(eta*eta)))
  endif else begin
     cross[i] = 0.
  endelse
endfor

return, cross
end

function kaulakys_crossh, A, nstarlow, llow, nupp, lupp, dE, E, cofm=cofm, version=version
; spin-averaged cross section for H collision on target

L_singlet = 5.965     ; scattering lengths in au, from Schwartz 1961 (agree with Bhatia 2007)
L_triplet = 1.7686
L_eff = sqrt(0.25 * (L_singlet^2.) + 0.75 * (L_triplet^2.))

if not keyword_set(version) then begin
   cross = kaulakys_cross(A, nstarlow, llow, nupp, lupp, 1.008, L_eff, dE, E, cofm=cofm)
endif else begin
   if version eq 1 then begin
   cross = kaulakys_cross(A, nstarlow, llow, nupp, lupp, 1.008, L_eff, dE, E, cofm=cofm)
   endif   
   if version eq 2 then begin
   cross = kaulakys_cross2(A, nstarlow, llow, nupp, lupp, 1.008, L_eff, dE, E, cofm=cofm)
   endif   
endelse

return, cross
end

function kaulakys_rateh, A, nstarlow, llow, nupp, lupp, dE, T, fast=fast
; the rate coefficient for collision with H
;
; where:
; A is the atomic mass of target (e.g. 40 for Ca)
; n and l are quantum numbers for upper and lower states of the target
;   nstar though, is the effective principal quantum number
; dE is the energy spacing of the levels in **eV**
; T is the temperature in K  -  may be an array
; result in cgs
;
;  fast keyword uses approximate method < sig(v) v > = sig(<v>) <v>
;  Kaulakys 1986 claims max error is 22%
;  

if keyword_set(fast) then begin
    kT = T *1.38066e-23/4.35981e-18 ;  kT in au  = T in au (k=1)
    mp = 1.008
    mu = mp*A/(A+mp) * 1822.89
    v = sqrt(8.*kT/!pi/mu) 
    E = 0.5*mu*v*v
    Q = kaulakys_crossh(A, nstarlow, llow, nupp, lupp, dE/27.2, E, cofm = 1) *(5.29177e-9^2)  ; convert to cm^2
    C = Q * v * 2.18769e8   ; convert v to cm/s
endif else begin
    ; this is brute force approach, perhaps would be better with adaptive method
    E = (dindgen(1000) * 0.05d0 + dE) / 27.2   ; 0->50eV in au
    Q = kaulakys_crossh(A, nstarlow, llow, nupp, lupp, dE/27.2, E, cofm = 1) 

    integ_rate, E*27.2, Q*(5.29177e-9^2), A, 1.008, T, C, /cofm
    ;plot, E*27.2, Q*(5.29177e-9^2), charsize = 3
endelse    
    
return, C
end