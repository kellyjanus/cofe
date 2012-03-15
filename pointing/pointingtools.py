import numpy as np

def NRAO_plot():
    ha=np.arange(0,12,.1)
    r=dict([(dec, np.degrees(compute_parallactic_angle((np.radians(ha/24.*360.)), np.radians(34), np.radians(dec)))) for dec in [-20, 0, 20, 30, 35, 40, 60,80]])
    rk = {}
    for k,v in r.iteritems():
        rk[k] = np.copy(v)
        rk[k][v<0]=v[v<0]+180
    from pylab import *
    figure()
    for dec,pa in rk.iteritems():
        plot(ha,pa,label=str(dec))
    grid()

def compute_parallactic_angle(ha, latitude, declination):
    """Inputs and output in radians"""
    return np.arctan2(
        np.sin(ha)/(
            np.tan(latitude)*np.cos(declination)-np.sin(declination)*np.cos(ha)
                   ), 1.)

def altaz2ha(elevation, azimuth, latitude):
    return np.arctan2( 
        -np.sin(azimuth)*np.cos(elevation), 
        -np.cos(azimuth)*np.sin(latitude)*np.cos(elevation)+
            np.sin(elevation)*np.cos(latitude)
                     )


if __name__ == '__main__':
    NRAO_plot()
