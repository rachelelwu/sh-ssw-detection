#!/bin/bash


#Loop from July (7) to December (12)
for mon in {7..12}; do
	echo "Running for month $mon..."
	        python dev_zeof_eof.py --lv 1000 --mon $mon
	done

	echo "All runs completed."

