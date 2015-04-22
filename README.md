# Pyxis-Issues
This is the initial pipeline to calibrate, image and extract sources observed via the KAT-7. Three files are attached: the first one is the python script, which contains the main procedures from calibration to source extraction. The second file is configuration file (.conf) that keeps the neccessary option set to imaging, creating directories, and others. The third one is a profile file (.profile) that contains the parameters obtained from Meqtrees to do the calibaration. In the head of this file, a name is given to call in the script.

The script is executed in the terminal using: 
$ pyxis -j8 OUTDIR=Test cal_ms
