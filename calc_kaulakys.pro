pro calc_kaulakys, infile, outdir
; this wrapper takes an input file and calculates data from it	

; the input file should be "unmerged" i.e. states are still divided into their components with different cores and thus cfps
; this must be so, since the spin redistribution requires the spin of the core to be known
;
; March 2023
; updated to deal with positive ions (previously neutrals only) - might break backward compatibility
; 

time = systime(1)

q = file_test(outdir)
if q then begin
    q1 = ' '
    print, outdir+' exists'
    print, ' w = wipe and redo '
    print, ' other = stop '
    read, q1, prompt = '(w/other): '
    case q1 of
        'w': begin
           print, ' Wiping ' + outdir
           spawn, 'rm -rf ' + outdir
        end
        else: begin
           print, ' Stopping.'
           stop
        end
    endcase
endif 
spawn, 'mkdir -p ' + outdir


Tmin = 50
Tmax = 1000
Tstep = 50

nT1=fix((Tmax-Tmin)/Tstep)+1
T1 = indgen(nT1)*Tstep + Tmin

Tmin = 500
Tmax = 20000
Tstep = 500

nT2=fix((Tmax-Tmin)/Tstep)+1
T2 = indgen(nT2)*Tstep + Tmin

T = [T1, T2]
nT = nT1 + nT2

T = T2 
nT = nT2

; read data
openr, lunm, infile, /get_lun
readf, lunm, ns, mass, ion, format = '(i10,f12.4,i10)'
label = strarr(ns)
g = intarr(ns)
E = dblarr(ns)
n = intarr(ns)
l = intarr(ns)
nstar = dblarr(ns)
s = intarr(ns)
sc = intarr(ns)
ionic = intarr(ns) * 0
core = strarr(ns)
cfp = fltarr(ns)
Hexc = intarr(ns) * 0

for i = 0, ns-1 do begin
    ;--
    ; ama51: original script to read input file
    ;--
	;ii=0
	;t1 = ' '
	;t2 = 0
	;t3 = 0.d0
	;t4 = 0
	;t5 = 0
	;t6 = 0.d0
	;t7 = 0.
	;t8 = 0.
	;t9 = 0
	;t10 = ' '
	;t11 = 0.
    ;readf, lunm, ii, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, format='(i3,a12,i3,f12.3,2i3,f12.3,2f5.1,i3,a10,f8.3)'

    ;--
    ; ama51: now use white space to separate inputs 
    ; (warning: do not use white space in your labels!)
    ;--
    str=''
    readf,lunm,str
    str=strsplit(str,/extract)
    label(i) = str[1]
    g(i) = double(str[2])
    E(i) = double(str[3])
    n(i) = double(str[4])
    l(i) = double(str[5])
    nstar(i) = double(str[6])
    s(i) = double(str[7])
    sc(i) = double(str[8])
    ionic(i) = double(str[9])
    core(i) = str[10]
    cfp(i) = double(str[11])
    ;--
    ; ama51: end modification
    ;--
    if strpos(label[i], 'H*') ge 0 then Hexc(i) = 1
endfor

close, lunm
free_lun, lunm

icov = ~ionic

iu = where(E ne shift(E, 1), nu)
E_uniq = E(iu)
indu = intarr(ns)
for i = 0, ns-1 do indu(i) = where(E_uniq eq E(i))

;stop

; write out the merged states list
   
openw, lunm, outdir+'/merged_states.txt', /get_lun
for i = 0, nu-1 do printf, lunm, i+1, label(iu(i)), E(iu(i)), format='(i3,a12,f12.3)'
close, lunm
free_lun, lunm


C_uniq = dblarr(nu,nu,nT)
Cdata = dblarr(ns,ns,nT)

; Kaulakys is for excitation - stated energy is transferred to the electron,  p_t is momentum threshold
; and so downward rates come from detailed balance.

ntrans = 0.5*(ns*ns - ns) 
dt_store = dblarr(ntrans)
trem_store = dblarr(ntrans)
ttot_store = dblarr(ntrans)
window, 2, xsize = 1000, ysize = 1000
!p.multi = [0,1,2]

count = 0
count_nonzero = 0
ttot = 0.
for i = 0, ns-2 do begin
   for j = i+1, ns-1 do begin
       print, ' '
       print, count+1, ' of ', fix(ntrans)
       print, i+1, j+1, nstar[i], nstar[j], l[i], l[j] 
       if sc[i] ne sc[j] then begin
         print, 'SKIPPING - different cores'
         goto, skip
       endif
       if core[i] ne core[j] then begin
         print, 'SKIPPING - different cores'
         goto, skip
       endif
       if ionic[i] or ionic[j] then begin
         print, 'SKIPPING - at least one state ionic'
         goto, skip
       endif
       if Hexc[i] or Hexc[j] then begin
         print, 'SKIPPING - excited H state involved'
         goto, skip
       endif       ; model can't couple different H excitations
       if (nstar[i]-l[i]) lt 0.1 then begin
         print, 'SKIPPING - Coulomb model breakdown at small n*-l'
         goto, skip
       endif  ; the Coulomb model breaks down
       if (nstar[j]-l[j]) lt 0.1 then begin
         print, 'SKIPPING - Coulomb model breakdown at small n*-l'
         goto, skip
       endif  ; the Coulomb model breaks down

       ; where multiple (two) parents but no equivalent electrons, we skip the second 
       ; to avoid counting twice.  This could in principle be accounted for via 
       ; angular momentum coupling to the core, see below

       if i ne 0 then begin
          if (label[i] eq label[i-1]) and (E[i] eq E[i-1]) and (cfp[i] eq 1) and (cfp[i-1] eq 1) then goto, skip
       endif 

       dE = (E[j]-E[i])/8065.45d0
       ts = {A:mass, ION:ion, N:n[i], L:l[i], ND:n[j], LD:l[j], NSTAR:nstar[i], DE:dE}
       ; since Emax = 30*kT in this routine, even with log grid, using << 30 pts is dangerous
       ;C = kaulakys_rateh(T, ts, method=2, npts=100, scat=1)
       ; for Eu II, 30 seems enough to get same result as 100 from 50 to 20000 K
       C = kaulakys_rateh(T, ts, method=2, npts=100, scat=0, delta_t = dt)
       ;C = kaulakys_rateh(T, ts, method=2, npts=100, scat=0, delta_t = dt)

       ; estimate total time remaining
       ttot = ttot + dt
       trem = (ttot / (count_nonzero+1.d0)) * ((count_nonzero+1.d0)/ (count+1.d0)) * (ntrans - count) / 86400.
       print, 'Estimated time remaining:', trem, ' days'
       
       ; plot time evolution
       ttot_store[count_nonzero] = ttot
       dt_store[count_nonzero] = dt 
       trem_store[count_nonzero] = trem
       plot, ttot_store[0:count_nonzero], dt_store[0:count_nonzero], xtitle = 'time[s]', ytitle = 'time [s]'
       plot, ttot_store[0:count_nonzero]/86400., trem_store[0:count_nonzero], xtitle = 'time[days]', ytitle = 'estimated time remaining [days]'

       	if sc[i] ne 0 then begin   ; in this case only one possible spin state

            ; find the other possible spin state
            sdash = sc[i] - 0.5    
            if sdash eq s[i] then sdash = sdash + 1

            ; calculate factor for spin change case
            fac = (2.*sdash+1.)/(2.*(2.*sc[i]+1.))*0.392

            ; apply factors appropriately
            if s[i] ne s[j] then begin  ; spin change
               C = fac*C
            endif else begin   ; no spin change
               C = (1-fac)*C
            endelse      

    	  endif	

    	C = C * cfp[i]^2.   ; since the cross section is proportional to the square of the initial state momentum-space wavefunction in Kaulakys (1985)
                          ; note that if sum(cfps^2)=1 as it should, then we will get back the original value below when we merge the components (different cores)
                          ; this also assumes the angular momentum coupling coeffs sum to 1 (sum over all projections).  This is something to prove in future.
                          ; for the moment one can at a minimum say angular momentum coupling is ignored

      ; an exception is cases where multiple (two) parents but no equivalent electrons
      ; see above
      
                          
                          

       ;facT = exp(dE/(8.617d-5*T))
       ;ind = where(finite(facT))
       ;Cdown = C * 0.d0
       ;Cdown[ind] = C[ind] * g[i]/g[j] * facT[ind]

       ; check for negative rate coefficients
       ; these have been seen in some cases with very small rate coeffs
       ; this is (hopefully) just when numerical precision is lost
       ; ideally one should do this at the cross section level, but it doesn't really matter.
       ind = where(C lt 0.d0, nind)
       if nind gt 0 then C[ind] = 0.d0

       Cdown = exp(alog(C *g[i]/g[j])+dE/(8.61733034d-5*T))   ; this limits numerical error at small T

       Cdata(i,j,*) = C
       Cdata(j,i,*) = Cdown

       print, 'T [K]'
       fmt = '(' + strcompress(string(nT, '(i6)'), /remove_all) + '(1X,i10))'
       fmt = '(10(1X,i10))'
       print, T, format = fmt
       fmt = '(' + strcompress(string(nT, '(i6)'), /remove_all) + '(1X,E10.2))'
       fmt = '(10(1X,E10.2))'
       print, 'rate coeff [cgs]'
       print, C, format = fmt
       print, 'rate coeff downwards [cgs]'
       print, Cdown, format = fmt
       count_nonzero= count_nonzero+1
       skip:
       count= count+1
   endfor
endfor

; merge all components

C_uniq = C_uniq * 0.d0
for i = 0, ns-1 do begin
	for j = 0, ns-1 do begin
        ii = indu(i)
        jj = indu(j)
        C_uniq(ii, jj, *) = C_uniq(ii, jj, *) + Cdata(i, j, *)
        if ii eq jj then C_uniq(ii, jj, *) = 0.d0
    endfor
endfor

; output


for iT = 0, nT-1 do begin
   sT = strcompress(string(T(iT), '(i6)'), /remove_all)
   fmt = '(' + strcompress(string(nu, '(i6)'), /remove_all) + '(1X,E13.3e3))'
   openw, lunm, outdir + '/' + sT + '_K.rates', /get_lun
   printf, lunm, C_uniq(*,*,iT), format = fmt
   close, lunm
   free_lun, lunm
endfor

print, 'input file: ', infile
print, 'output directory: ', outdir
print, 'number of transitions: ', count
print, 'number of non-zero transitions: ', count_nonzero
print, 'time for calc_kaulakys: ', systime(1) - time
print, 'time per transition: ', (systime(1) - time)/(count_nonzero*1.)


end

