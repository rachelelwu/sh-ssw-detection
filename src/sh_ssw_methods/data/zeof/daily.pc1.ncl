; USAGE:> ncl lv=50 mon=6 daily.pc1.ncl

begin

 eday = (/31,29,31,30,31,30,31,31,30,31,30,31/)

 print(eday(mon-1))

 hlat = -90
 llat = -20

 samplef = addfile("daily.z1.197906.nc","r")
 sampled = short2flt(samplef->z(0,:,:)) 
 printVarSummary(sampled)

 clat = sqrt(cos(sampled&latitude*3.14/180))
 printVarSummary(clat)

 sampleDim = dimsizes(sampled)
 print(sampleDim)

 fileNames = systemfunc("ls daily.z"+lv+".{1979,198*,199*,200*,201*,202[0-3]}"+sprinti("%2.2i",mon)+".nc")

 inFile = addfiles(fileNames,"r")
 ListSetType(inFile,"join")
 z = inFile[:]->z

 printVarSummary(z)
 
 gphI = z/9.8
 copy_VarCoords(z,gphI)

 printVarSummary(gphI) 

 delete(z)

 gphI!0 = "year"
 gphI&year = ispan(1979,2023,1)
 gphI!1 = "day"
 gphI!2 = "lat"
 gphI!3 = "lon"

 anom = gphI - conform(gphI,dim_avg_n_Wrap(gphI,0),(/1,2,3/))

 copy_VarCoords(gphI,anom)
 printVarSummary(anom) ; [year | 45] x [day | 30] x [lat | 71] x [lon | 360]

; print(anom&lat+"  "+clat)

 pwt = anom*conform(anom,clat,2)
 pwt@_FillValue = -999

 copy_VarCoords(anom,pwt)
 printVarSummary(pwt)

 ncdf = addfile ("ERA5.EOF1.z"+lv+".Mon"+sprinti("%2.2i",mon)+".nc","r")
 ev = ncdf->eof1
; 
 printVarSummary(ev)
 
 ts = anom(:,:,0,0)

 printVarSummary(ts)  ; [year | 45] x [day | 30]

 yr = ispan(1979,2023,1)

 nyr = dimsizes(gphI&year)

 do y=0,nyr-1
  do d=0,eday(mon-1)-1
   ts(y,d) = sum(pwt(y,d,:,:)*ev(:,:))
  end do
 end do

 ts@_FillValue = -999
  
 sts = dim_standardize_n_Wrap(ts,1,0)

 stddv = dim_stddev_n_Wrap(ts,0)

 print(stddv(24))
 print(ts(:,0)+"  "+sts(:,0))

 strdata = new(nyr,string)
 strdata = yr

 do y=0,nyr-1
  do d=0,eday(mon-1)-1
    strdata(y) = strdata(y) + sprintf("%9.3f",sts(y,d))
  end do
  asciiwrite("daily.pc1.z"+lv+"."+yr(y)+sprinti("%2.2i",mon)+".txt",strdata(y))
 end do
end
