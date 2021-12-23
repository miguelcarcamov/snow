residual_ms=$1
model_fits=$2
file_restored=$3
file_residuals=$4
weighting=$5
robust=$6

echo $residual_ms $model_fits $file_restored $weighting $robust
/home/simon/anaconda3/bin/python pyrarestore.py $residual_ms $model_fits $file_restored $file_residuals $weighting $robust

echo back to shell

# source ~/gitcommon/objectoriented_selfcal/casa6/bin/activate
# python --version
