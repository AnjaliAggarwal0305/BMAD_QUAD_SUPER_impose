#!/bin/bash
#rm -r ../production
mk clean
mk
../production/bin/test_orbit_correction 
#python COSY_closed_orbit.py
#python x,y_vs_time.py
#python sx_sy_sz_vs_time.py
