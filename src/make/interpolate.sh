
###
#
# Interpolate paths to binaries from bins.cfg into target scripts
# R references need to be placed in: 
# lib/analyze2.R
# lib/libyogiplot.R
# lib/master.R
# lib/slave.R
# lib/synergizer.R
#

config=$1
buildDir=$2

#read bin locations from cfg.bin
declare -A bins=(
	["Rbin"]=`grep Rscript $config|cut -d, -f2` 
	["BLASTbin"]=`grep blast $config|cut -d, -f2` 
	["SGEbin"]=`grep Rscript $config|cut -d, -f2`
)

function interpolate {
	
	targetfile=$1

	#Check if file exists
	if [ ! -r $targetfile ] 
	then
		echo "Cannot read file $targetfile!"
		exit 1
	fi

	#Load file contents into variable
	contents=`cat $targetfile`

	#iterator over remaining arguments
	for key in ${!bins[@]}
	do
		#escape slash characters (don't even think about messing with this!)
		value=`echo ${bins[$key]}|sed -e 's/\//\\\\\//g'`
		#replace current variable with its value
		contents=`echo "$contents"|sed "s/\\\$$key/$value/g"`
		# echo "$key => ${bins[$key]}"
	done

	#return results
	echo "$contents">$targetfile
}

interpolate "$buildDir/lib/analyze2.R"
interpolate "$buildDir/lib/master.R"
interpolate "$buildDir/lib/slave.R"
