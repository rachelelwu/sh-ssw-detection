import numpy as np
import pandas as pd

from scipy.ndimage import uniform_filter1d

def drophf4n(
    hf,
    SDY1,
    EDY2,
    DTL,
    DUL,
    NYRS,
    SY,
    NDYS,
    NDSEP,
    NEGHF=None
):
    """
    Line-by-line translation of IDL DROPHF4N

    ; HEATFLUX criterion
    ; tjr 12/1/21
    """

    # equivalent of IDL THRESH, TIME, NDSEPA
    THRESH = DUL
    TIME = DTL
    NDSEPA = NDSEP
    
    # if n_elements(NEGHF) eq 0 then NEGHF=0
    if NEGHF is None:
        NEGHF = 0

    hf1 = hf.copy()  # local copy
 
    T = DTL

    # DAT=[[ ],[ ]]  in IDL means 2x0 empty # output array(4,n): 4 entries (year, day, 0, HFsum) per event, n events
    # easiest in python: list of lists
    DAT = []

    NDW = EDY2 - SDY1 + 1 # number of days in winter
    if NDW <= 0:
        NDW += NDYS

    # nan = fltarr(NDYS)
    nan = np.zeros(NDYS, dtype=np.float32)  # 1 year of 0s

    ML = 80  # number of days to search after a 0 crossing for SSW criteria (vortex return)

    if EDY2 < 180: # this must be NH 
        # hfu = [nan, hf1]
        hfu = np.concatenate([nan, hf1]) # unsmoothed hf
        SY1 = SY - 1
        NYRS1 = NYRS + 1
    else:
        hfu = hf1
        SY1 = SY
        NYRS1 = NYRS

    # for YR = 0l, NYRS1-2 do begin ; tjr 4/23/20
    for YR in range(0, NYRS1 - 1):
        
        mx = min(NYRS1 * NDYS - 1, (YR + 2) * NDYS + SDY1) # 2 full years or less

        hyu = hfu[(YR * NDYS + SDY1):(mx + 1)]
        # print(hyu)

        # W = fltarr(NDW*2)
        W = np.zeros(NDW * 2, dtype=np.float32) # weights
        for ii in range(0, NDW * 2):
            W[ii] = np.exp(-1.0 * ii / TIME)

        i = 0 # days since SDY1
        j = 0 # day of last event

        # repeat begin ... endrep until
        while i <= NDW - 1:
            # print(i)
            
            # hfsum = total(reverse(hyu(j:i)) * W(0:i-j),/nan)
            subarray = hyu[j:i+1][::-1]
            weights = W[0:i - j + 1]
            hfsum = np.nansum(subarray * weights) # weighted sum from prior event (j) to day i

            if NEGHF:
                hfsum *= -1

            if hfsum >= THRESH: # event found
                if NEGHF: # negative event
                    while i < NDW - 1 and hyu[i] < 0:
                        i += 1
                else: # positive event
                    while i < NDW - 1 and hyu[i] > 0:
                        i += 1
                i -= 1 # i: last day of hf > 0

                # weighted sum of prior hf
                subarray = hyu[j:i+1][::-1]
                weights = W[0:i - j + 1]
                hfsum = np.nansum(subarray * weights)

                # in some rare cases it happens that hfsum < THRESH by advancing i
                if NEGHF:
                    if -hfsum < THRESH:
                        hfsum = -THRESH
                else:
                    if hfsum < THRESH:
                        hfsum = THRESH

                k = i - 3 # arbitrary adjustment so that SLP crosses 0 at onset
                if j != 0 and (k - j) < NDSEPA: # ignore event if previous event was less than 20 days ago
                    break  # jump out of repeat

                y1 = YR
                d1 = k + SDY1
                if d1 >= NDYS and SDY1 > EDY2:
                    d1 -= NDYS
                    y1 += 1
                # append like [[DAT],[new]]
                DAT.append([y1 + SY1, d1, 0, hfsum])

                j = i
            i += 1

    # in IDL it returns DAT, so
    return np.array(DAT).T  # same shape as IDL


def read_daily_hf(FN):
    """
    # ; read daily heat flux file in binary 
    # ; in the following, what I loosely call 'heatflux' are actually vertical EP-fluxes at 100 hPa

    Input(s):
    FN (str) - filepath to heat flux binary file from ERA5 - 'f2bym_100_20_90_SH_ERA5_1979_2023.bin'
    NDYS (int) - number of day of years in the data
    NYRS (int) - number of years in the data
    NM (int) - number of modes in the data

    Output(s):
    hfo (array) - heat flux from binary file
    """
    
    with open(FN, "rb") as f:
        NDYS, NYRS, NM = np.fromfile(f, dtype=np.int32, count=3)
        hfo = np.fromfile(f, dtype=np.float32, count=NDYS*NYRS*NM)

    return hfo, NDYS, NYRS, NM

def detrend_hf_anom(hfo, NDYS, NYRS, NM):
    """
    # ; detrended heat flux anomalies

    Input(s):
        hfo (array) - daily heat flux data (in dim: NDYS*NYRS*NM)
        NDYS (int) - number of day of years in the hf data
        NYRS (int) - number of years in the hf data
        NM (int) - number of modes in the hf data
    
    Output(s):
        hf_mode0 (array) - detrended hf anomalies, HF(16425), mode 0 only (= unfiltered sum of all Foruier modes)

    """
    from scipy.ndimage import uniform_filter1d
    
    # 1. scaling
    HF = np.reshape(hfo, (NDYS, NYRS, NM), order="F") / -1e5
    
    # 2. daily climatology (mean over years)
    HFclim = np.nanmean(HF, axis=1)  # shape (NDYS, NM)

    # 3. smooth climatology twice
    for m in range(NM):
        HFclim[:, m] = uniform_filter1d(HFclim[:, m], size=25, mode="wrap")
        HFclim[:, m] = uniform_filter1d(HFclim[:, m], size=25, mode="wrap") #  smooth twice, works quite nicely

    # 4. daily anomalies
    for i in range(NYRS):
        HF[:,i,:] -= HFclim

    # 5. remove slowly varying trend (not sure how important the next two steps are)
    for m in range(NM):
        for d in range(NDYS):
            window = HF[max(d-10,0):min(d+10+1,NDYS),:,m]
            ss = np.nanmean(window, axis=0)
            s = uniform_filter1d(ss, size=20, mode="nearest")
            HF[d,:,m] -= s

    # 6. remove linear trend
    b = np.arange(NYRS)
    
    for m in range(NM):
        for d in range(NDYS):
            y = HF[d,:,m]
            coeffs = np.polyfit(b, y, 1)
            trend = np.polyval(coeffs, b)
            HF[d,:,m] -= trend
    
    # 7. flatten mode 0
    hf_mode0 = HF[:,:,0].reshape((NDYS*NYRS), order="F") # HF(16425), mode 0 only (= unfiltered sum of all Foruier modes)

    return hf_mode0


def main_EPflux_reichler(FN, DT=None, DU=None, output_csv=True):

    
    # 0. inputs
    SY = 1979
    EY = 2023

    MODEL = 'ERA5'
    
    NYRS = EY - SY + 1
    
    
    # detection criteria: dates (1/1 = 0)
    NDSEP = 20			#; separation in days between individual events: Charlton et al. 2007
    LATITUDE = -60			#; SH
    HLEVEL = 100 			#; FZ level in hPa (so far only for NNR and ERA5)
    NDYS = 365
    SDY1 = 120; M0 = 4; D0 = 0	#; start day for detection, 5/1	corresponds to 11/1 over NH
    EDY2 = NDYS-1			#; end day for detection (last day of year)

    if DT is None:
        DT = 60
    if DU is None:
        DU = 36
    
    # 1. read heat flux from binary file
    hfo, NDYS, NYRS, NM = read_daily_hf(FN)

    # 2. detrend heat flux anomalies
    hf_mode0 = detrend_hf_anom(hfo, NDYS, NYRS, NM)

    # 3. detect events
    DAT = drophf4n(hf_mode0, SDY1, EDY2, DT, DU, NYRS, SY, NDYS, NDSEP, NEGHF=0)

    if output_csv:
        file = f"FZ{HLEVEL}_SH_{MODEL}_tau{DT}_sfz{DU}.csv"

        with open(file, "w") as fout:
            fout.write(f"# Data: {MODEL} {SY}-{EY}\n")
            fout.write(f"# Definition: accumulated weighted (tau = {DT} d) FZ{HLEVEL} > {DU} WDU\n")
            fout.write(f"# Search between day {SDY1} and day {EDY2}.\n")
            fout.write(f"# Events are the same if occurring within {NDSEP} days.\n")
            fout.write(f"# Reference: Reichler and Jucker 2022 (WCD)\n")
        
            fout.write("year,month,day\n")
        
            ndpm = [31,28,31,30,31,30,31,31,30,31,30,31]  # number of days per month
            dates = []

        
            if DAT.shape[1] == 0:
                print("No events to write.")
            else:
                for i in range(DAT.shape[1]):
                    doy = DAT[1, i]
                    year = int(DAT[0, i])
        
                    # convert day-of-year to month/day
                    k = 0
                    while doy >= ndpm[k]:
                        doy -= ndpm[k]
                        k += 1
                        if k >= 12:
                            break  # safety check
        
                    month = k + 1  # IDL month numbering starts at 0
                    day = doy + 1
                    fout.write(f"{year},{month},{day}\n")

                    # build datetime64
                    year = int(DAT[0, i])
                    doy  = int(DAT[1, i])   # assuming DOY is integer-like
                     # leap-year safe: Jan 1 + doy days
                    d = np.datetime64(f"{year:04d}-01-01") + np.timedelta64(doy, "D")
                    dates.append(d)

        
        print(f"{file} saved.")

    dates = np.array(dates, dtype="datetime64[D]")  # daily precision


    return dates, DAT