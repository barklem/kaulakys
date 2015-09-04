n = 9
; energies of states in /cm
ELi = [0., 14903., 27206., 30925., 31283., 35012., 36469., 36623., 36630.]
SLi = ['2s',   '2p',    '3s',   '3p',   '3d',   '4s',   '4p',   '4d',   '4f']
nLi = [2,2,3,3,3,4,4,4,4]
lLi = [0,1,0,1,2,0,1,2,3]

T  = (dindgen(9)+1)*1000.

openw, 1, 'kaul_Li.txt'
printf, 1, n, n_elements(T)
printf, 1, T, format = '(9f10.0)'
flush, 1


for i = 0, n-3 do begin
   for j = i+1, n-2 do begin
       print, i, j
       nstarlow = sqrt(109678./(41449.-ELi(i)))
       ts = {A:6.941, N:nLi[i], L:lLi[i], ND:nLi[j], LD:lLi[j], NSTAR:nstarlow, DE:(ELi[j]-ELi[i])/8065.45d0}
       C = kaulakys_rateh(T, ts, method=2, npts=10)
       printf, 1, SLi[i],  ' ',SLi[j]
       printf, 1, C, format = '(9e10.3)'
       flush, 1
   endfor
endfor
close, 1

end