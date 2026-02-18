; USAGE:> ncl lv=50 event_selection.ncl

begin

 system("rm -f out.pc1.z"+lv+".*txt")

 thrshld = 1
if(lv.eq.50) then
 duration = 9      ; total duration = first detected day + duration     
else
 duration = 13
end if

 eday = (/30,31,31,30,31,30,31/)  ; from June to December

 sday = new(7,integer)
 accday = new(7,integer)

 do m=0,6
  accday(m) = sum(eday(0:m))
  sday(m) = accday(m) - eday(m)
 end do
 print(sday+"  "+accday)

 tot_day = accday(6)

 year = ispan(1979,2023,1)
 nyr = dimsizes(year)

 data = new((/nyr,accday(6)/),float) 

 printVarSummary(data)

;----------------------------------
; these do loops are to generate an array of year x day of pc1 values from the files written by years and months
; 
 do y=0,nyr-1
  do m=0,6
   mon = m+6
   f = asciiread("daily.pc1.z"+lv+"."+year(y)+sprinti("%2.2i",mon)+".txt",(/1, eday(m)+1/),"float")
   data(y,sday(m):accday(m)-1) = f(0,1:eday(m))
   delete(f)
  end do
 end do

 data@_FillValue = -999
 print(data(0,:))

;--------------------------------
; to set date information
;
 days = ispan(1,tot_day,1)
 stime1 = new((/nyr,tot_day/),integer)
 stime_yy = new((/nyr,tot_day/),integer)
 stime_mm = new((/nyr,tot_day/),integer)
 stime_dd = new((/nyr,tot_day/),integer)

 month = days
 month = 6

 hour = days
 hour = 0

 minute = days
 minute = 0

 second = days
 second = 0

 units = "days since 1979-01-01 00:00:00"

 do y=0,nyr-1
  yyy = days
  yyy = year(y)    ; randomly chosen leap year

  date = ut_inv_calendar(yyy,month,days,hour,minute,second,units,0)

  stime = ut_calendar(date,0)
  printVarSummary(stime)

  stime1(y,:) = year(y)*10000 + floattoint(stime(:,1)*100) + floattoint(stime(:,2))
  stime_yy(y,:) = year(y)
  stime_mm(y,:) = floattoint(stime(:,1))
  stime_dd(y,:) = floattoint(stime(:,2))
 end do

; print(stime1(0,:)+"  "+stime_yy(0,:)+"  "+stime_mm(0,:)+"  "+stime_dd(0,:))

; time coordinate values are set
;----------------------------------------------------------
; Identifying positive ZEOF1 events above the threshold
;
 count = new((/nyr,tot_day-duration/),integer)
 count = 0

 count2 = new((/nyr/),integer)
 count2 = 0

 do y=0,nyr-1
  do d=0,(tot_day-1-duration)

   if(data(y,d).gt.thrshld) then
    subIDX = ind(data(y,d+1:d+duration).gt.thrshld) 
    if(dimsizes(subIDX).eq.duration) then
;      print("first date of an event "+stime1(y,d)+"  "+data(y,d))
      count(y,d) = 1
    end if
    delete(subIDX)
   end if

  end do
 end do

 do y=0,nyr-1

   xidx = ind(count(y,:).eq.1)
   print(xidx)

   nxidx = dimsizes(xidx)
;   print(year(y)+"  "+nxidx)   

   if(nxidx.eq.1.and.ismissing(xidx)) then
     print("missing "+xidx)
     delete([/xidx,nxidx/])
   else

;-----------------------------------------------
;  write the first event detected
;-----------------------------------------------
     asciiwrite("out.pc1.z"+lv+"."+stime1(y,xidx(0))+".txt",stime_yy(y,xidx(0))+","+stime_mm(y,xidx(0))+","+stime_dd(y,xidx(0)))

;----------------------------------------------
; Find independent events separated by 60 days and write them out
;
     count2(y) = 1

     do j=0,nxidx-2
      del = xidx(j+1)-xidx(j)
       if(del.ge.60) then
;       print(xidx(j)+"  "+xidx(j+1))
        asciiwrite("out.pc1.z"+lv+"."+stime1(y,xidx(j+1))+".txt",stime_yy(y,xidx(j+1))+","+stime_mm(y,xidx(j+1))+","+stime_dd(y,xidx(j+1)))
        count2(y) = count2(y) + 1
       end if
     end do
     delete([/xidx,nxidx/])

   end if

 end do

 print("Total number identified "+count2)

;---------------------------------
; concatenate the out.aao.zLV*.txt to a file to have the final EOF1.ZLV.gt.1.txt
;---------------------------------
 system("cat out.pc1.z"+lv+".*txt > EOF1.Z"+lv+".gt.1.txt") 
 system("rm -f out.pc1.z"+lv+".*txt")

end
