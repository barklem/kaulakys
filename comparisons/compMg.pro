lcao = read_ascii('MgH_LCAO_6000_K_cut.rates', data_start=1)
quant = read_ascii('MgHquantum6000', data_start=1)
kaul = read_ascii('MgH_kaulakys_scat=1_6000_K_cut.rates', data_start=1)

lcao_data = lcao.field1
quant_data = quant.field1
kaul_data = kaul.field1

n=8

; Extract the lower and upper triangular arrays.
l = fltarr(n, n)
u = fltarr(n, n)
FOR j = 1,n - 1 DO l[0:j-1,j] = lcao_data[0:j-1,j]
FOR j=0,n - 1 DO u[j:*,j] = lcao_data[j:*,j]
lcao_down = u
lcao_up = l

l = fltarr(n, n)
u = fltarr(n, n)
FOR j = 1,n - 1 DO l[0:j-1,j] = quant_data[0:j-1,j]
FOR j=0,n - 1 DO u[j:*,j] = quant_data[j:*,j]
quant_down = u
quant_up = l

l = fltarr(n, n)
u = fltarr(n, n)
FOR j = 1,n - 1 DO l[0:j-1,j] = kaul_data[0:j-1,j]
FOR j=0,n - 1 DO u[j:*,j] = kaul_data[j:*,j]
kaul_down = u
kaul_up = l

; cut ionic

lcao_down = lcao_down[0:n-2, 0:n-2]
lcao_up = lcao_up[0:n-2, 0:n-2]
quant_down = quant_down[0:n-2, 0:n-2]
quant_up = quant_up[0:n-2, 0:n-2]
kaul_down = kaul_down[0:n-2, 0:n-2]
kaul_up = kaul_up[0:n-2, 0:n-2]

psymcircle
plot, quant_up, lcao_up / quant_up, psym=8, /ylog, /xlog
oplot,  quant_up, (lcao_up + kaul_up) / quant_up, psym = 4


compratio = lcao_up / quant_up
indcomp = where(finite(compratio) eq 1 and compratio gt 0.d0, ncomp)
meanratio = mean(compratio(indcomp)) ; mean ratio
stddevratio = stddev(compratio(indcomp))  ; stddev
weights = quant_up(indcomp)
sumweights = total(weights)
wmeanratio = total(compratio(indcomp)*weights)/sumweights
wstddevratio=sqrt(total((compratio(indcomp)-wmeanratio)^2*weights)/sumweights*ncomp/(ncomp-1))

print, 'lcao/quant'
print, meanratio, stddevratio
print, wmeanratio, wstddevratio



compratio = (lcao_up + kaul_up) / quant_up
indcomp = where(finite(compratio) eq 1 and compratio gt 0.d0, ncomp)
meanratio = mean(compratio(indcomp)) ; mean ratio
stddevratio = stddev(compratio(indcomp))  ; stddev
weights = quant_up(indcomp)
sumweights = total(weights)
wmeanratio = total(compratio(indcomp)*weights)/sumweights
wstddevratio=sqrt(total((compratio(indcomp)-wmeanratio)^2*weights)/sumweights*ncomp/(ncomp-1))


print, '(lcao + kaul)/quant'
print, meanratio, stddevratio
print, wmeanratio, wstddevratio


plot, quant_up[0:2,0:2], lcao_up[0:2,0:2] / quant_up[0:2,0:2], psym=8, /ylog, /xlog
oplot,  quant_up[0:2,0:2], (lcao_up[0:2,0:2] + kaul_up[0:2,0:2]) / quant_up[0:2,0:2], psym = 4

plot, quant_up[0:2,0:2]+ kaul_up[0:2,0:2], lcao_up[0:2,0:2] / quant_up[0:2,0:2], psym=8, /ylog, /xlog
oplot,  quant_up[0:2,0:2]+ kaul_up[0:2,0:2], (lcao_up[0:2,0:2] + kaul_up[0:2,0:2]) / quant_up[0:2,0:2], psym = 4

stop

end