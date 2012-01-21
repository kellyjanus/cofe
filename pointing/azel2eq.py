import pycfitsio as fits

freq = 15
pointing_file = fits.open('/home/zonca/COFE/data/Level1/0.5/v5_all_%dghz_pointing.fits' % freq)
data_file = fits.open('/home/zonca/COFE/data/Level1/0.5/all_%dGHz_v0.5_data.fits' % freq)
servo_file = fits.open('/home/zonca/COFE/data/Level1/0.5/utservo_all_v0.5.fits')

# pointing channel is the column in the pointing file
pnt_ch = 1

MISSIONSTART = 16.8 #from altitude
MISSIONEND = 36.84 #from issue with latitude

ut = data_file['TIME'].read_column('UT') 
good = (ut > MISSIONSTART) & (ut < MISSIONEND)

az = pointing_file[0].read_column('az%d' % pnt_ch)[good]
el = pointing_file[0].read_column('el%d' % pnt_ch)[good]

servo_ut = servo_file['GYRO_HID'].read_column('UT')
servo_good = (servo_ut > MISSIONSTART) & (servo_ut < MISSIONEND)
servo_alt = servo_file['GYRO_HID'].read_column('HYBRIDALTITUDE')
good_alt = servo_good & (servo_alt > 4000) & (servo_alt < 1e5)
servo_lat = servo_file['GYRO_HID'].read_column('HYBRIDLATITUDE')
good_lat =  servo_good &(servo_lat < .65) & (servo_lat > .56)
servo_lon = servo_file['GYRO_HID'].read_column('HYBRIDLONGITUDE')
good_lon = servo_good & (servo_lon > -2) & (servo_lon < -1.5)
