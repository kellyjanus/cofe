import matplotlib.pyplot as plt
import sys
sys.path.append('../code')

from lib import *

base_folder ='/COFE'
day = 'all'
version='debug'

gyro = fits.read(os.path.join(base_folder, 'servo', '%s.fits' % day), 'GYRO_HID')

plt.figure(); plt.title('Gyro computer clock')
plt.plot(gyro['COMPUTERCLOCK'], label='before offsets')

offsets = find_clock_offsets_from_gpstime(gyro['COMPUTERCLOCK'], gyro['GPSTIME'])

plt.plot(gyro['COMPUTERCLOCK'], label='after offsets')
plt.legend(loc=0); plt.xlabel('Sample index')

utcc, ut = create_ut(gyro)

plt.figure(); plt.title('UT')
plt.plot(ut)
plt.legend(loc=0); plt.xlabel('Sample index')

plt.figure(); plt.title('UT vs CC')
plt.plot(utcc, ut)
plt.legend(loc=0); plt.xlabel('computerClock')

servo_file = os.path.join(base_folder,'servo','%s.fits' % day)
create_utservo(servo_file, offsets, utcc, ut)

plt.figure(); plt.title('UT of devices')
for device in DEVICES:
    plt.plot(fits.read('utservo.fits', device)['UT'], label=device)

plt.legend(loc=0); plt.xlabel('Sample index')

for freq in [10, 15]:

    revcounter = fits.read(os.path.join(base_folder, 'servo', '%s.fits' % day), REVCOUNTER_LABEL[freq])
    sci_file = os.path.join(base_folder,str(freq),'%s.fits'%day)
    data = fits.read(sci_file)

    plt.figure()
    plt.title('%d GHz REV' % freq)
    plt.plot(data['REV'],'.',label='SCI REV before fix')
    plt.plot(revcounter['VALUE'] ,'.',label='SERVO REV before fix')

    servo_range = remove_reset(revcounter['VALUE'], offsetsci=data['REV'][0])

    # apply offsets to revcounter cc
    revcounter['COMPUTERCLOCK'] = apply_cc_offsets(revcounter['COMPUTERCLOCK'], offsets)

    # oversample revcounter cc and value to 140 Hz in order to interpolate over gaps
    uniform_rev_cc = np.arange(revcounter['COMPUTERCLOCK'][servo_range][0], revcounter['COMPUTERCLOCK'][servo_range][-1], (1/140.)*2e9, dtype=np.double)

    plt.plot(revcounter['VALUE'][servo_range] ,'.',label='SERVO REV after fix')

    uniform_rev = np.interp( uniform_rev_cc, revcounter['COMPUTERCLOCK'][servo_range].astype(np.double), revcounter['VALUE'][servo_range].astype(np.double))
    plt.plot(uniform_rev ,'.',label='SERVO OVERSAMPLED')
    plt.legend(loc=0); plt.xlabel('Sample index')

    sci_cc = create_utscience(sci_file, gyro, revcounter, offsets, utcc, ut, freq)
    plt.figure()
    plt.title('%d GHz sci CC' % freq)
    plt.plot(sci_cc)
    plt.legend(loc=0); plt.xlabel('Sample index')

    create_sync_servo(servo_file, offsets, utcc, ut, sci_cc, freq)

    #plt.figure()
    #plt.title('%d GHz MAG' % freq)
    #plt.plot(cctotime(sci_cc), np.degrees(np.arctan2(splitted_data['CHANNEL_15']['T'], splitted_data['CHANNEL_14']['T'])),label='SCI')
    #plt.plot(cctotime(synched_data['TIME']['COMPUTERCLOCK']), np.degrees(np.arctan2(synched_data['ANALOGCHANNELS']['CHANNEL_31'], synched_data['ANALOGCHANNELS']['CHANNEL_30'])),label='ANALOGCH')
    #plt.legend(); plt.xlabel('hours')
