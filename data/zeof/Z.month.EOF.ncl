;*************************************************
; Usage e.g.: ncl  mon=6 lv=1 Z.month.EOF.ncl
;************************************************
begin

  print(lv+"  "+mon)

  longMonthNames = new((/12/),string)
  longMonthNames(:) = (/"January","February","March","April","May","June","July","August","September","October","November","December"/)

  fileName = systemfunc("ls ../../GPH/z"+lv+".era5.{1979,198[0-9],199[0-9],200[0-9],201[0-9],202[0-3]}??.nc")
  print(fileName)

  inFile = addfiles( fileName, "r" )
  ListSetType( inFile, "cat" )
  z = lonFlip(inFile[:]->z(:,{-90:-20:2},::2))  

  z@_FillValue = -999

  printVarSummary(z)

  print(z&latitude)

  z!1 = "lat"
  z!2 = "lon"

  OBmslp_1d = z
  printVarSummary(OBmslp_1d)  ; [time | 540] x [lat | 71] x [lon | 360]

;  print(OBmslp_1d&lat)

  wgt = sqrt(cos(OBmslp_1d&lat*0.01745329))
  printVarSummary(wgt)

  pwt = OBmslp_1d*conform(OBmslp_1d,wgt,1)/9.8
  pwt@_FillValue = -999
  printVarSummary(pwt)

  copy_VarCoords(OBmslp_1d,pwt)
  printVarSummary(pwt)

  pwt1mon = pwt((mon-1)::12,:,:)
  printVarSummary(pwt1mon)

  clim = dim_avg_n_Wrap(pwt1mon,0)

  printVarSummary(clim)
  print(max(clim)+"  "+min(clim))

  pwtanom = pwt1mon - conform(pwt1mon,clim,(/1,2/))

  printVarSummary(pwtanom)

  delete(z)
  delete(OBmslp_1d)
  delete(pwt)

  printVarSummary(pwt1mon)

  copy_VarCoords(pwt1mon,pwtanom)

  delete(pwt1mon)

  pwtanom!0 = "time"

  pwtanom&time = ispan(1979,2023,1)

  pwtanom@_FillValue = -999

  printVarSummary(pwtanom)

;;------------------------------------------------------
;  EOF analysis

  neval = 2

  nlat = dimsizes(pwtanom&lat(:))
  nlon = dimsizes(pwtanom&lon)

  ev = new((/neval,nlat,nlon/),float)
  printVarSummary(ev)

  ev = eofunc(pwtanom(lat|:,lon|:,time|:),neval,False) 
  printVarSummary(ev)

  print(" % variance  "+ev@pcvar)
  print("EVAL  "+ev@eval)

  ev!0 = "eigenv"
  ev!1 = pwtanom!1
  ev!2 = pwtanom!2

  ev&eigenv = ispan(1,neval,1)
  ev&lat = pwtanom&lat
  ev&lon = pwtanom&lon
;
  printVarSummary(ev)

  nyr = dimsizes(pwtanom&time)
  ev_ts = eofunc_ts_Wrap(pwtanom(lat|:,lon|:,time|:),ev,True)
;
  print(ev_ts(0,:))
  printVarSummary(ev_ts)
 
  if(ev(0,{-60},{90}).lt.0) then
    n_ev_ts = dim_standardize(ev_ts(0,:),1)*-1.0
    eof1 = ev(0,:,:)*-1.0
  else
    n_ev_ts = dim_standardize(ev_ts(0,:),1)
    eof1 = ev(0,:,:)
  end if

  copy_VarCoords(ev(0,:,:),eof1)

  EOFp = regCoef(n_ev_ts,pwtanom(lat|:,lon|:,time|:))

  copy_VarCoords(pwtanom(0,:,:),EOFp)

  system("/bin/rm -f ERA5.EOF1.z"+lv+".Mon"+sprinti("%2.2i",mon)+".nc")
  ncdf = addfile ("ERA5.EOF1.z"+lv+".Mon"+sprinti("%2.2i",mon)+".nc","c")
  ncdf->eof1 = eof1

;  asciiwrite("Z."+lv+"hPa.PC1.1979-2023.Mon"+sprinti("%2.2i",mon)+".txt",n_ev_ts)
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("x11" ,"Z."+lv+"hPa.EOF1.Mon"+sprinti("%2.2i",mon))        

  gsn_define_colormap(wks,"BlueWhiteOrangeRed")        ; choose colormap

  res                       = True     ; plot mods desired

  res@gsnPolar = "SH"

  res@gsnDraw = False 
  res@gsnFrame = False
  res@vpWidthF = 0.7

  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn off contour lines

  res@gsnContourZeroLineThicknessF = 3
  res@gsnContourNegLineDashPattern = 1

  res@cnLineLabelsOn = False
  res@cnInfoLabelOn         = False
  res@lbLabelBarOn          = True
  res@lbLabelFontHeightF  = 0.015      ; make labels smaller
  res@lbLabelStride         = 2        ; label interval

  res@lbTitleString = "m"
  res@lbTitlePosition = "Right"
  res@lbTitleDirection = "Across"

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF        = -175.
  res@cnMaxLevelValF        = 175.
  res@cnLevelSpacingF       = 25.0     ; contour spacing
  res@gsnSpreadColors       = True     ; use full range of color map

  res@pmTickMarkDisplayMode = "Always" ; use NCL default lat/lon labels

  res@gsnAddCyclic          = True

  res@mpCenterLonF  = 0
  res@mpMinLatF            = -90
  res@mpMaxLatF            = -20.
  res@mpGridLatSpacingF = 30
  res@mpGridLonSpacingF = 45
  res@gsnPolarLabelFontHeightF = .011
  res@mpGeophysicalLineColor = "gray50"

  res@gsnCenterString = longMonthNames(mon-1)
  res@gsnCenterStringFontHeightF = 0.025
  res@gsnRightString = "EOF 1: "+sprintf("%4.1f",ev@pcvar(0))+"%"
  res@gsnRightStringFontHeightF = 0.02
  res@gsnLeftString = lv+"hPa"
  res@gsnLeftStringFontHeightF = 0.02

  res@gsnLeftStringFontHeightF = 0.02

  plot = gsn_csm_contour_map_polar(wks,EOFp(:,:), res)

  res_c = True
  res_c@cnLinesOn = True
  res_c@cnFillOn = False

  res_c@gsnDraw = False
  res_c@gsnFrame = False

;  res_c@cnMaxLevelValF = 50000.
;  res_c@cnMinLevelValF = -500.
;  res_c@cnLevelSpacingF = 100.

  res_c@gsnLeftString = " "
  res_c@gsnRightString = " "
  res_c@gsnContourLineThicknessF = 2.
  res_c@gsnContourZeroLineThicknessF = 3
  res_c@gsnContourNegLineDashPattern = 2
  res_c@gsnAddCyclic          = True

  plot1 = gsn_csm_contour(wks,clim,res_c)

  overlay(plot,plot1)

  draw(plot)
  frame(wks)

end
