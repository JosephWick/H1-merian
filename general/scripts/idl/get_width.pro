; FUNCTION GET_WIDTH
;
; Returns the width of the HI line at some percentage of the peak.
;
; INPUTS:
;
;    SPECTRUM: The spectrum of the HI (ie from total_line_profile)
;
;    VAXIS: The velocity axis corresponding to the specturm (ie from
;           total_line_profile
;
;    WIDTH: The percentage of the peak at which to find the width.
;
;
; OUTPUTS:
;
;    W: The width at chosen percentage of the peak.
;
;
function get_width, spectrum, vaxis, width, v1, v2, widthval, vnet=vnet, onepeak=onepeak

if not (keyword_set(spectrum)) then begin
    print,'get_width(spectrum, [vaxis, width, v1, v2, widthval])'
    return,-1
endif


; check inputs.
if (not keyword_set(spectrum)) or (not keyword_set(vaxis))$
  or not keyword_set(width) then begin
    print,'Usage: getWidth(spectrum, velocity_axis, width)'
    return,-1
endif
if (n_elements(spectrum) ne n_elements(vaxis)) then begin
    print,'Spectrum must have the same number of elements ' + $
      'as the velocity axis'
    return,-1
endif
if (width ge 1) then begin
    print,"WIDTH must be less than 1 (percentage of peak)."
    return,-1
endif

if not keyword_set(vnet) then vnet = 0.

; adjust velocities just in case.
;vave = total(vaxis * spectrum) / total(spectrum)
;vaxisp = vaxis - vave

peak = max(spectrum)
widthval = peak * width
i1 = where(vaxis lt 0)

; figure out first side
; make sure you're looking past the peak
indlow = max(spectrum[where(vaxis le vnet)], indmax)
blah = min(abs(spectrum[0:indmax] - widthval), i1)
;blah = min(abs(spectrum[where(vaxis le vnet)] - widthval), i1)
if ((spectrum[i1] gt widthval) and (i1 gt 0))then begin
    m = (spectrum[i1] - spectrum[i1-1]) / (vaxis[i1] - vaxis[i1-1])
    b = spectrum[i1] - m * vaxis[i1]
    v1 = (widthval - b) / m
endif else if ((spectrum[i1] lt widthval) and $
               (i1 lt n_elements(vaxis)-1)) then begin
    m = (spectrum[i1+1] - spectrum[i1]) / (vaxis[i1+1] - vaxis[i1])
    b = spectrum[i1] - m * vaxis[i1]
    v1 = (widthval - b) / m
endif else v1 = vaxis[i1]

; figure out second side
; make sure you're looking past the peak
indhi = max(spectrum[where(vaxis gt vnet)], indmax)
blah = min(abs(spectrum[min(where(vaxis gt vnet))+indmax:n_elements(spectrum)-1] - widthval), i2)
i2 = i2 + min(where(vaxis gt vnet))+indmax
if ((spectrum[i2] gt widthval) and (i2 gt 0)) then begin
    m = (spectrum[i2] - spectrum[i2-1]) / (vaxis[i2] - vaxis[i2-1])
    m = (spectrum[i2] - spectrum[i2-1]) / (vaxis[i2] - vaxis[i2-1])
    b = spectrum[i2] - m * vaxis[i2]
    v2 = (widthval - b) / m
endif else if ((spectrum[i2] lt widthval) and $
               (i2 lt n_elements(vaxis)-1)) then begin
    m = (spectrum[i2+1] - spectrum[i2]) / (vaxis[i2+1] - vaxis[i2])
    b = spectrum[i2] - m * vaxis[i2]
    v2 = (widthval - b) / m
endif else v2 = vaxis[i2]
;    endelse
w = abs(v2 - v1)
print,"Width at " + strtrim(ceil(100*width),1) + "% of the peak is " + $
  strtrim(w,1) + " km/s."



return,w

end
