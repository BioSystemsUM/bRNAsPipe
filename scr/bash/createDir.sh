# This script has a function that receives an array of directory paths
# and checks whether they exist, if not it creates them
createDirIfMissing()
{
	arr=($@)
	for dir in ${arr[@]}
	do [[ -d ${dir} ]] || mkdir ${dir}
	done
}
