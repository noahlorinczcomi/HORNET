git clone https://github.com/noahlorinczcomi/HORNET.git && cd HORNET
curl -o hornet_data.tar http://hal.case.edu/~njl96/hornet_data.tar.gz
tar -xf hornet_data.tar && del hornet_data.tar
mkdir tempfiles results plots
pip install -r requirements.txt
del example_man_vol_plots.png
msg /w %username% HORNET is now installed in the %cd% folder

