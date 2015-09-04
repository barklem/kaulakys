n = 9
; energies of states in /cm
ENa = [0.000, 16956., 25739., 29172., 30266., 33200., 34548., 34588., 35050.]; / 8065.45d0
SNa = ['3s',   '3p',    '4s',   '3d',   '4p',   '5s',   '4d',   '4f',   '5p']
nNa = [3,3,4,3,4,5,4,4,5]
lNa = [0,1,0,2,1,0,2,3,1]

T  = (dindgen(9)+1)*1000.

openw, 1, 'kaul_Na.txt'
printf, 1, n, n_elements(T) 
printf, 1, T, format = '(9f10.0)'
flush, 1

for i = 0, n-3 do begin
   for j = i+1, n-2 do begin
       print, i, j
       nstarlow = sqrt(109678./(41449.-ENa(i)))
       ts = {A:23., N:nNa[i], L:lNa[i], ND:nNa[j], LD:lNa[j], NSTAR:nstarlow, DE:(ENa[j]-ENa[i])/8065.45d0}
       C = kaulakys_rateh(T, ts, method=2, npts=10)
       printf, 1, SNa[i],  ' ',SNa[j]
       printf, 1, C, format = '(9e17.7)'
       flush, 1
   endfor
endfor
close, 1

end