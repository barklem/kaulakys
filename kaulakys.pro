; routines for calculation of Kaulakys free electron model
;
; Paul Barklem
; written Jan-March 2010
; lots of new parts Jan 2013
;
; Dec 2022 - March 2023 added extension to ions:
; ion is the total charge on the atom (0 = neutral, 1=singly ionised)
;
; Notes: 1. various parameters in this code may need tweaking for your particular case (esp. pmax and npts in Inl)
;
; 

pro gwfsq_gen, ion, nstar, l
; pregenerates the momentum space wavefunction

common WAVEFUNC, parr, gsqarr

npts = 10000  
pm   = 2.d0/!pi/nstar   ; <p>  n*=1 -> 0.64  n*=100 -> 6.4e-3
pmin = 1e-6
pmax = 1e2

; generate log grid for npts-1 points on log grid (more points at small momenta)

pstep = (alog(pmax) - alog(pmin))/(npts-2)
parr = exp(alog(pmin) + pstep*dindgen(npts-1))

; add zero
parr = [0.d0, parr]

gsqarr = gwfsq(ion, nstar, l, parr)

end

function gwfsq_int, ion, nstar, l, p, noint=noint

common WAVEFUNC, parr, gsqarr

if keyword_set(noint) then begin
  gsq = gwfsq(ion, nstar, l, p)
endif else begin  
  gsq = interpol(gsqarr, parr, p)
endelse

return, gsq
end  


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

function ikern1, x, private
p = private.p
x1 = private.x1
x2 = private.x2
ans = x*0.
ind = where((1-x)*(x2-x)*(x-x1) gt 0., nind)
if nind gt 0 then ans[ind] = scatamp_h(p,x[ind]) / sqrt((1-x[ind])*(x2-x[ind])*(x[ind]-x1))
ind2 = where(~finite(ans), nind2)
if nind2 gt 0 then stop
return, ans
end

function ikern2, x, private
p = private.p
x1 = private.x1
x2 = private.x2
ans = x*0.
ind = where((1-x)*(-x1-x)*(x2+x) gt 0., nind)
if nind gt 0 then ans[ind] = scatamp_h(p,x[ind]) / sqrt((1-x[ind])*(-x1-x[ind])*(x2+x[ind]))
ind2 = where(~finite(ans), nind2)
if nind2 gt 0 then stop
return, ans
end

function igr1, p, x1, x2, xt2
private = {x1:x1, x2:x2, p:p}
ans = qpint1d('ikern1', x1, xt2, private, breakpoints = [1., x1, x2], epsrel=1e-3)
return, ans
end

function igr2, p, x1, x2, xt1
private = {x1:x1, x2:x2, p:p}
ans = qpint1d('ikern2', -x2, xt1, private, breakpoints = [1., -x1, -x2], epsrel=1e-3)
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

function kkernel, p, private
pt = private.pt
x1 = private.x1
x2 = private.x2
nstar = private.nstar
l = private.l
ion = private.ion
return, ksigma(p, pt, x1, x2) * gwfsq_int(ion, nstar,l,p, noint=1) * p*p
end

function kaulakys_crossh_best, E, private, cofm=cofm

A = private.nstar
ion = private.ion
nstar = private.nstar
n = private.n
l = private.l
nd = private.nd
ld = private.ld
dE = private.dE

APert = 1.008
mu = A*APert/(A+APert)*1822.89    ; reduced mass for target-pertuber in au
if keyword_set(cofm) then begin
   Ecm = E
endif else begin
   Ecm = double(A)/(A+APert) * E
endelse

; allow array of energies
cross = dblarr(n_elements(Ecm))
for i = 0, n_elements(Ecm)-1 do begin
  v = sqrt(2.d0*Ecm[i]/mu)                         ; relative velocity in au
  pt = abs(dE) / (2.d0*v)                          ; momentum transfer in au 
  if Ecm[i] gt dE then begin                       ; note dE < 0 for exothermic
     
   x1 = (sqrt(n*n-(l+0.5d0)^2.)*sqrt(nd*nd-(ld+0.5)^2.)-(l+0.5)*(ld+0.5))/(n*nd)
   x2 = (sqrt(n*n-(l+0.5d0)^2.)*sqrt(nd*nd-(ld+0.5)^2.)+(l+0.5)*(ld+0.5))/(n*nd)
   
   inf = !values.d_infinity
   private2 = {pt:pt, x1:x1, x2:x2, nstar:nstar, l:l, ion:ion}
   igral = qpint1d('kkernel', pt, +inf, private2, epsrel=1e-3)
   cross[i] = (2.*ld+1.)*sqrt(2.)/2./nd^5./v^2. * igral
   
   
   if (~finite(cross[i])) then begin
   npts= 1000 & pmin = pt & pmax = 2.
   p = dindgen(npts) * (pmax-pmin)/(npts-1) + pmin
   ii = kkernel(p, private2)
   stop
   endif
   
  endif else begin
     cross[i] = 0.d0
  endelse
endfor

return, cross
end

function Inl_kernel, p
; calculates the kernel of the Inl integral 

common QUANTUM, nstarp, lp, ion
common SHARE, pt

wfsq = gwfsq_int(ion,nstarp,lp,p,noint=1)
kernel = wfsq*(p-pt)*p

return, kernel
end

function Inl, ion, nstar, l, pt
; calculates the Inl function, a type of overlap integral

common QUANTUM, nstarp, lp, ionp
nstarp = nstar
lp = l
ionp = ion

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

function kaulakys_cross, A, ion, nstar, l, nd, ld, APert, Lscat, dE, E, cofm=cofm 
; the cross section in scattering length approximation (eq. 18, 1991 paper)
; where:
; A is the atomic mass of target (e.g. 40 for Ca)
; ion is the total charge on the atom (0 = neutral, 1=singly ionised)
; n and l are quantum numbers for upper and lower states of the target
; nstar though, is the effective principal quantum number
; APert is the atomic mass of the perturber (e.g. 1 for H)
; Lscat is the scattering length of the perturber-electron system in au
; dE is the energy spacing of the levels in au
; E is the kinetic energy in the laboratory (target's) frame in au
;   if cofm is set then it's in the centre-of-mass frame

common SHARE, pt

mu = A*APert*1822.89/(A+APert)                   ; reduced mass for target-pertuber
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
     cross[i] = 2.*!pi*Lscat^2.*(2.*ld+1.)/nd^5./v^2. * Inl(ion, nstar, l, pt)
  endif else begin
     cross[i] = 0.
  endelse
endfor

return, cross
end

function kaulakys_cross2, A, ion, nstar, l, nd, ld, APert, Lscat, dE, E, cofm=cofm 
; the cross section in scattering length approximation (eq. 12, 1991 paper)
; assume degeneracy
; where:
; A is the atomic mass of target (e.g. 40 for Ca)
; ion is the total charge on the atom (0 = neutral, 1=singly ionised)
; n and l are quantum numbers for upper and lower states of the target
; nstar though, is the effective principal quantum number
; APert is the atomic mass of the perturber (e.g. 1 for H)
; Lscat is the scattering length of the perturber-electron system in au
; dE is the energy spacing of the levels in au
; E is the kinetic energy in the laboratory (target's) frame in au
;   if cofm is set then it's in the centre-of-mass frame

common SHARE, pt

mu = A*APert*1822.89/(A+APert)                   ; reduced mass for target-pertuber
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
     Igral =  Inl(ion, nstar, l, 0.)
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

function kaulakys_cross86, A, ion, nstarlow, llow, nupp, APert, Lscat, dE, E, cofm=cofm 
; this is for nl -> n', from 1986 paper
; using the Inl calculated directly here

common SHARE, pt

mu = A*APert*1822.89/(A+APert)                   ; reduced mass for target-pertuber
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
     cross[i] = 2.*!pi*Lscat^2./nupp^3./v^2. * Inl(ion, nstarlow, llow, pt)
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

mu = A*APert*1822.89/(A+APert)                   ; reduced mass for target-pertuber
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

function kaulakys_crossh, E, private, cofm=cofm, version=version
; spin-averaged cross section for H collision on target

A = private.A
ion = private.ion
nstarlow = private.nstar
nlow = private.n
llow = private.l
nupp = private.nd
lupp = private.ld
dE = private.dE

L_singlet = 5.965     ; scattering lengths in au, from Schwartz 1961 (agree with Bhatia 2007)
L_triplet = 1.7686
L_eff = sqrt(0.25 * (L_singlet^2.) + 0.75 * (L_triplet^2.))

if not keyword_set(version) then begin
   cross = kaulakys_cross(A, ion, nstarlow, llow, nupp, lupp, 1.008, L_eff, dE, E, cofm=cofm)
endif else begin
   if version eq 1 then begin
   cross = kaulakys_cross(A, ion, nstarlow, llow, nupp, lupp, 1.008, L_eff, dE, E, cofm=cofm)
   endif   
   if version eq 2 then begin
   cross = kaulakys_cross2(A, ion, nstarlow, llow, nupp, lupp, 1.008, L_eff, dE, E, cofm=cofm)
   endif   
endelse

return, cross
end


function ratekernel_best, E, private, kT=kT
return, kaulakys_crossh_best(E, private, cofm = 1) * E * exp(-E/kT) 
end

function ratekernel, E, private, kT=kT
return, kaulakys_crossh(E, private, cofm = 1) * E * exp(-E/kT) 
end

function kaulakys_rateh, T, tstruct, method=method, npts=npts, plt=plt, scat=scat, delta_t=delta_t
; the rate coefficient for collision with H
;
; T is the temperature in K  -  may be an array
;
; tstruct is a structure detailing the transition and should include:
; A - the atomic mass of target (e.g. 40 for Ca)
; ion is the total charge on the atom (0 = neutral, 1=singly ionised)
; n and l - quantum numbers for initial state of the target
; nd and ld - quantum numbers for final state of the target
; nstar - the effective principal quantum number of the **initial** state
; dE is the energy spacing of the levels in **eV**
;    sign is important : >0 endothermic, <0 exothermic
;
; result in cgs
;
;  method keyword :
;  1           : uses fast approximate method < sig(v) v > = sig(<v>) <v>
;                Kaulakys 1986 claims max error is 22%
;  2           : use brute force on prescribed log E grid with npts points 
;                default is 10 points
;                plt turns on plotting for checking grid
;  3 (default) : use adaptive integration routine 
;                this will be very slow for multiple T, since adaptive part depends on T
;                and so cross sections can't be reused
;
;   scat -> use scattering length approximation

A      = tstruct.A
ion    = tstruct.ion
n      = tstruct.n
l      = tstruct.l
nd     = tstruct.nd
ld     = tstruct.ld
nstar  = tstruct.nstar
dE     = tstruct.dE / 27.2116   ; -> au

private = {A:A, ion:ion, nstar:nstar, n:n, l:l, nd:nd, ld:ld, dE:dE}

if ~keyword_set(method) then method = 3
if ~keyword_set(npts) then npts = 10

time = systime(1)

kT = T *1.38066e-23/4.35981e-18 ;  kT in au  = T in au (k=1)

; pregenerate initial state wavefunction
;gwfsq_gen, nstar, l

case method of
   1: begin
         mp = 1.008
         mu = mp*A/(A+mp) * 1822.89  
         v = sqrt(8.*kT/!pi/mu) 
         E = 0.5*mu*v*v
         if keyword_set(scat) then begin
            Q = kaulakys_crossh(E, private, cofm = 1) *(5.29177e-9^2)  ; convert to cm^2
         endif else begin
            Q = kaulakys_crossh_best(E, private, cofm = 1) *(5.29177e-9^2)  ; convert to cm^2
         endelse   
         C = Q * v * 2.18769e8   ; convert v to cm/s
      end 
   2: begin
         Emin = max([1d-6, dE])                     ; in au
         Emax = max([kT * 30., dE * 30.])
         Estep = (alog(Emax) - alog(Emin))/(npts-1)
         E = exp(alog(Emin) + Estep*dindgen(npts))  ; grid in au
         E = [E, Emax*10, Emax*100]                 ; extend grid roughly at large E
         if keyword_set(scat) then begin
            Q = kaulakys_crossh(E, private, cofm = 1) 
         endif else begin
            Q = kaulakys_crossh_best(E, private, cofm = 1) 
         endelse
         if keyword_set(plt) then begin
            plot, E*27.2116, Q, /xlog
            oplot, [kT, kT]*27.2116 , [1e-30, 1e10], linestyle = 2
            oplot, [Emax, Emax]*27.2116 , [1e-30, 1e10], linestyle = 2
         endif
         integ_rate, E*27.2116, Q*(5.29177e-9^2), A, 1.008, T, C, /cofm
      end
   3: begin
         mp = 1.008
         mu = mp*A/(A+mp) * 1822.89  
         Emin = max([0.d0, dE])
         C1 = sqrt(8/!dpi/mu) / (kT^1.5) 
         inf = !values.d_infinity
         nT = n_elements(T)
         igral = dblarr(nT)
         for i = 0, nT-1 do begin
            if keyword_set(scat) then begin
              igral[i] = qpint1d('ratekernel', Emin, +inf, private, functargs = {kT:kT[i]}, epsrel=1e-3)
            endif else begin
              igral[i] = qpint1d('ratekernel_best', Emin, +inf, private, functargs = {kT:kT[i]}, epsrel=1e-3)
            endelse  
         endfor
         C = C1 * igral * (5.29177e-9^2) * 2.18769e8   ; au -> cgs
      end
endcase
delta_t = systime(1) - time
print, 'time for kaulakys_rateh: ', delta_t
return, C
end