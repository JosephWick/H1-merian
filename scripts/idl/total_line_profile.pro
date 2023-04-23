; PRO TOTAL_LINE_PROFILE
;
; Creates a line profile for all the detected emission. Also finds W20
; if requested. Will hanning smooth if requested.
;
; INPUTS:
;
;    CUBE: the data array as read in by mrdfits (read_cube_fits.pro)
;
;    HEADER: the header structure as read in by read_cube_fits.pro
;
;    VAXIS: will be set to the velocity axis.
;
;    SPECTRUM: will be set to the spectrum/intensity axis (in msun
;    at least)
;
;
; OPTIONAL INPUTS:
;
;    HANNING: integer describing the type of hanning smoothing to do.
;
;
;  
pro total_line_profile, cube, header, vaxis, spectrum, w, $
                        hanning=hanning, boxcar=boxcar, $
                        width=width, doplot=doplot, outfile=outfile, $
                        name=name, interactive=interactive, click=click,$
                        onepeak=onepeak, shift=shift, bin=bin


if not (keyword_set(cube) and keyword_set(header)) then begin
    print,"total_line_profile, cube, header, [vaxis, spectrum,"
    print,"                    hanning=hanning, boxcar=boxcar,"
    print,"                    width = width, /plot, outfileoutfile]"
    return
endif
; make sure input is okay.
if not keyword_set(cube) then begin
    print,"CUBE must be set to a data cube array."
    return
endif else if n_elements(size(cube,/dimensions)) ne 3 then begin
    print,"CUBE must be a three-dimensional array."
endif
if not keyword_set(header) then begin
    print,"HEADER keyword not set; spectrum will be in arbitrary units."
endif
if (keyword_set(hanning) and keyword_set(boxcar)) then begin
    print,"set either HANNING or BOXCAR for smoothing, but not both."
    return
endif
if keyword_set(hanning) then begin
    if ((hanning mod 2) eq 0) then begin
        print,"Smoothing must be an odd number."
        return
    endif
endif
if keyword_set(boxcar) then begin
    if ((boxcar mod 2) eq 0) then begin
        print,"Smoothing must be an odd number."
        return
    endif
endif
if not keyword_set(width) then width=0.2
if (width ge 1) then begin
    print,"WIDTH must be less than 1 (percentage of peak)."
    return
endif




; set up velocity axis and scale the cube using the header.
if keyword_set(header) then begin
    vaxis = dindgen(header.naxis3) * header.cdelt3 + header.crval3
endif else begin
    vaxis = dindgen(n_elements(spectrum))
endelse


; pretend we're a radio telescope..
spectrum = total( total(cube, 1, /double), 1, /double)


; bin if requested
; bin keyword is binsize in km/s
if keyword_set(bin) then begin

    ; figure out what the new size of the array will be
    newsize = n_elements(spectrum) * header.cdelt3 / bin

    ; resize the spectrum and the velocity axis
    spectrum = congrid(spectrum, newsize)
    vaxis = congrid(vaxis, newsize)

endif


; smooth if requested
if keyword_set(hanning) then begin
    kernel = fltarr(hanning) + 1.
    for i=0,floor(hanning/2.) do begin
        kernel[i] = 2.^i
        kernel[hanning - 1 - i] = 2.^i
    endfor
    kernel = kernel / double(total(kernel))
    spectrum = convol(spectrum, kernel,/edge_truncate)
endif else if keyword_set(boxcar) then begin
    kernel = (fltarr(hanning) + 1.) / double(hanning)
    spectrum = convol(spectrum, kernel, /edge_truncate)
endif


; adjust velocities
if keyword_set(shift) then begin
    vave = total(vaxis * spectrum) / total(spectrum)
    vaxis = vaxis - vave
endif


; match peaks:
if keyword_set(onepeak) then begin
;    blah = max(spectrum, i)
;    vnet = vaxis[i]
    vnet = total(spectrum * vaxis) / total(spectrum)
    v1 = vnet - 100
    v2 = vnet + 100
endif else begin
    blah = max(spectrum[where(vaxis le 0)], i1)
    v1 = vaxis[i1]
    blah = max(spectrum[where(vaxis gt 0)], i2)
    i2 = i2 + n_elements(where(vaxis le 0))
    v2 = vaxis[i2]
    vnet = (v2 + v1)/2
endelse
if keyword_set(shift) then vaxis = vaxis - vnet

if keyword_set(width) then begin
    w = get_width(spectrum, vaxis, width, v1w, v2w, widthval, vnet=vnet, onepeak=onepeak)
endif

if keyword_set(outfile) or keyword_set(doplot) then begin
    if not keyword_set(doplot) then begin
        set_plot,'ps'
        device,filename = outfile + '.lp.ps'
    endif

    ; set the title of the graph
    titl = 'HI Line Profile'
    if keyword_set(name) then titl = titl + ', ' + name

    ; figure out axis range
    xmin = 100.*ceil(v1/100.)-100
    xmax = 100.*ceil(v2/100.)+100
    xrang = [-1. * max(abs([xmin,xmax])), max(abs([xmin,xmax]))]
;    xrang=[-1000,1000]
openw, lun, 'hi.data.dat', /get_lun
for p = 0L,n_elements(spectrum)-1 do printf, lun, vaxis[p], spectrum[p]
close, lun
free_lun, lun
    plot, vaxis, spectrum, xtitle='Velocity (km/s)', pos=[.2,.2,.9,.9], $
      ytitle='M' + textoidl('_{HI}') +' (M' + textoidl('_{sun}') + ')',$
      title=titl, $
      xrange = xrang, xstyle=1
    if keyword_set(width) then begin
        widthval = width * max(spectrum)
        oplot, [v1w, v2w], [widthval,widthval], psym=2
        oplot, [v1w, v2w], [widthval, widthval], linestyle=2

        xyouts, 0.15, 0.9, /normal, 'W' + strtrim(round(100*width),2) $
          +' / 2 = '+ strtrim(string(w/2., format='(d10.1)'),2)+ ' km/s', charsize=1.5
    endif
    
    if not keyword_set(doplot) then begin
        device,/close
        set_plot,'x'
    endif 
endif


end 
