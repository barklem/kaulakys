pro calc_kaulakys, infile, outdir
; this wrapper takes an input file and calculates data from it	


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

Tmin = 1000
Tmax = 20000
Tstep = 1000

nT=fix((Tmax-Tmin)/Tstep)+1
T = indgen(nT)*Tstep + Tmin


; read data
openr, lunm, infile, /get_lun
readf, lunm, ns, mass, format = '(i10,f12.4)'
label = strarr(ns)
g = intarr(ns)
E = dblarr(ns)
n = intarr(ns)
l = intarr(ns)
nstar = dblarr(ns)
s = intarr(ns)
sc = intarr(ns)
ionic = intarr(ns)
core = strarr(ns)
cfp = fltarr(ns)
   
for i = 0, ns-1 do begin
	ii=0
	t1 = ' '
	t2 = 0
	t3 = 0.d0
	t4 = 0
	t5 = 0
	t6 = 0.d0
	t7 = 0.
	t8 = 0.
	t9 = 0
	t10 = ' '
	t11 = 0.
    readf, lunm, ii, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, format='(i3,a12,i3,f12.3,2i3,f12.3,2f5.1,i3,a10,f8.3)'
    label(i) = t1
    g(i) = t2
    E(i) = t3
    n(i) = t4   
    l(i) = t5
    nstar(i) = t6
    s(i) = t7
    sc(i) = t8
    ionic(i) = t9
    core(i) = t10
    cfp(i) = t11
endfor

close, lunm
free_lun, lunm

icov = ~ionic

iu = where(E ne shift(E, 1), nu)
E_uniq = E(iu)
indu = intarr(ns)
for i = 0, ns-1 do indu(i) = where(E_uniq eq E(i))

; write out the merged states list
   
openw, lunm, outdir+'/merged_states.txt', /get_lun
for i = 0, nu-1 do printf, lunm, i+1, label(iu(i)), E(iu(i)), format='(i3,a12,f12.3)'
close, lunm
free_lun, lunm


C_uniq = dblarr(nu,nu,nT)
Cdata = dblarr(ns,ns,nT)

for i = 0, ns-2 do begin
   for j = i+1, ns-1 do begin
       print, i+1, j+1
       if sc[i] ne sc[j] then goto, skip
       if core[i] ne core[j] then goto, skip
       if ionic[i] or ionic[j] then goto, skip  
       if (nstar[i]-l[i]) lt 0.1 then goto, skip   ; the Coulomb model breaks down
       if (nstar[j]-l[j]) lt 0.1 then goto, skip

       dE = (E[j]-E[i])/8065.45d0
       ts = {A:mass, N:n[i], L:l[i], ND:n[j], LD:l[j], NSTAR:nstar[i], DE:dE}
       C = kaulakys_rateh(T, ts, method=2, npts=10, scat=1)
       ;C = kaulakys_rateh(T, ts, method=1)
       
       	if sc[i] ne 0 then begin

   			fac = (2.*s[j]+1.)/(2.*(2.*sc[i]+1.))*0.392
            if s[i] ne s[j] then begin
               C = fac*C
            endif else begin
               C = (1-fac)*C
            endelse      

    	endif	

    	C = C * cfp[i]^2.

       Cdown = C * g[i]/g[j] * exp(dE/(8.617e-5*T))
       Cdata(i,j,*) = C
       Cdata(j,i,*) = Cdown

       fmt = '(' + strcompress(string(nT, '(i6)'), /remove_all) + '(1X,i10))'
       print, T, format = fmt

       fmt = '(' + strcompress(string(nT, '(i6)'), /remove_all) + '(1X,E10.2))'
       print, C, format = fmt
       print, Cdown, format = fmt
       skip:
   endfor
endfor

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

end

