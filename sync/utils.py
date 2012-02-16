import numpy as np

def cctotime(c):
    return (c-c[0])*2.e-9

def cctout(cc, gpstime):
    """cc and gpstime must be already synched"""
    utc = np.mod(((gpstime+15.)/3600.), 24)
    ut = utc[0]+(((cc - cc[0])*2e-9)/3600.)
    return ut

def fix_counter_jumps(d):
    """Removes jumps by linear fitting and removing points further than a predefined threshold, then linearly interpolates to the full array"""
    THRESHOLD = 3
    t = np.arange(len(d))
    ar, br = np.polyfit(t,d,1)
    lin_d = np.polyval([ar, br], t)
    valid = np.abs(d - lin_d) < THRESHOLD
    d_fix = np.interp(t, t[valid], d[valid], left=d[valid][0]-1, right=d[valid][-1]+1)
    return d_fix.astype(np.int)

def fix_counter_jumps_diff(d):
    """Removes 1 sample jumps by checking the sample to sample difference"""
    THRESHOLD = 5
    t = np.arange(len(d))
    fix = np.abs(np.diff(d))<THRESHOLD
    lin_d = d[fix]
    d_fix = np.interp(t, t[fix], lin_d, left=lin_d[0]-1, right=lin_d[-1]+1)
    return d_fix.astype(np.int)

def remove_reset(d, offsetsci=None):
    """Gets longest time period between resets
    
    Parameters
    ----------
    d : ndarray
        revcounter array
    offsetsci : num
        start sample of science data

    Returns
    -------
    s : slice
        slice of good data range between resets
    """
    # first check for 20bit jumps
    jump20bit_indices, = np.where(np.logical_and(np.diff(d) < -2**20*.9, np.diff(d) > -2**20*1.1))
    print('20bit jumps at:')
    print(jump20bit_indices)
    for j in jump20bit_indices:
        d[j+1:] += 2**20
    if not offsetsci is None:
        d += 2**20 * np.round((offsetsci - d[0])/2**20)
    print('remove dip')
    start,=np.where(np.diff(d) < -500000)
    stop,=np.where(np.diff(d) > 500000)
    if len(start)>0:
        d[start[1]-1:stop[1]+1] = -1
    reset_indices, = np.where(np.diff(d) < -300000)
    real_reset = []
    for i in reset_indices:
        if not (1.5e6<i<1.6e6): 
        #if abs(d[i+2]-d[i]) >= abs(d[i+1]-d[i]) and not (1.5e6<i<1.6e6): 
            #it is a single point, not a reset
            real_reset.append(i)
    reset_indices = np.array(real_reset)

    print('reset jumps at:')
    print(reset_indices)
    sections_boundaries = np.concatenate([[0], reset_indices +1 ,[len(d)]])
    sections_lengths = np.diff(sections_boundaries)
    max_len = sections_lengths.argmax()
    s = np.zeros(len(d), dtype=np.bool)
    start = np.max([150000, sections_boundaries[max_len]])
    stop = np.min([sections_boundaries[max_len+1], d.searchsorted(4865400)])
    s[start+1:stop-1] = True
    s[d==-1] = False
    s[d==0] = False
    single_sample_jump, = np.where(np.diff(d[:stop-1])<0)
    assert np.all(np.diff(single_sample_jump)>2)
    s[single_sample_jump + 1] = False
    print('selecting section with no resets from index %d to index %d %.2f %% of the total length' % (start, stop, (stop-start)*100./len(d)))
    return s

def index_of_single_sample_jumps(d):
    jumps,=np.where(np.abs(np.diff(d))>.005 )
    single_jumps = jumps[np.diff(jumps)==1] + 1
    return single_jumps

def make_monotonic(d):
    jumps, = np.where(np.diff(d)<0)
    print('Fixing jumps at indices %s, at relative position %s' % (str(jumps), str(jumps.astype(np.double)/len(d)) ))
    for j in jumps:
        d[j+1:] += d[j] - d[j+1]
