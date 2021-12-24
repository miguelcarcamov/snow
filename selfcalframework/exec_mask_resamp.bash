file_canvas=$1
file_mask=$2
file_out=$3

temp=$( realpath "$0"  )
zdir=$(dirname "$temp")

/home/simon/anaconda3/bin/python $zdir/mask_resamp.py $file_canvas $file_mask $file_out

echo back to shell

# source ~/gitcommon/objectoriented_selfcal/casa6/bin/activate
# python --version
