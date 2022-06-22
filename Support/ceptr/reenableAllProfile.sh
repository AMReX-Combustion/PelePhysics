declare -a fileArray=("ceptr/symbolic_math.py"\
                      "ceptr/jacobian.py")

length=${#fileArray[@]}

# Iterate the string array using for loop
for (( i=0; i<${length}; i++ ));
do
    file="${fileArray[$i]}"
    sed -i.bu 's/# @profile/@profile/g' $file
    rm $file.bu
done

