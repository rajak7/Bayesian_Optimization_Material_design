#!/bin/sh
i=0
submitted=0
for M1 in Mo W
do
  for X1 in S Se Te
  do
    for M2 in Mo W
    do
      for X2 in S Se Te
      do
        for M3 in Mo W
        do
          for X3 in S Se Te
          do
            symm_index=0
            echo $M1${X1}2-$M2${X2}2-$M3${X3}2
	    third=${M3}${X3}
	    second=${M2}${X2}
	    first=${M1}${X1}
	    if [[ $first == MoS ]]
	    then
		symm_index=$((symm_index+0))
	    elif [[ $first == MoSe ]]
            then
                symm_index=$((symm_index+1))
            elif [[ $first == MoTe ]]
            then
                symm_index=$((symm_index+2))
            elif [[ $first == WS ]]
            then
                symm_index=$((symm_index+3))
            elif [[ $first == WSe ]]
            then
                symm_index=$((symm_index+4))
	    else
		symm_index=$((symm_index+5))
	    fi
            if [[ $second == MoS ]]
            then
                symm_index=$((symm_index+0))
            elif [[ $second == MoSe ]]
            then
                symm_index=$((symm_index+6))
            elif [[ $second == MoTe ]]
            then
                symm_index=$((symm_index+12))
            elif [[ $second == WS ]]
            then
                symm_index=$((symm_index+18))
            elif [[ $second == WSe ]]
            then
                symm_index=$((symm_index+24))
            else
                symm_index=$((symm_index+30))
            fi
            if [[ $third == MoS ]]
            then
                symm_index=$((symm_index+0))
            elif [[ $third == MoSe ]]
            then
                symm_index=$((symm_index+36))
            elif [[ $third == MoTe ]]
            then
                symm_index=$((symm_index+72))
            elif [[ $third == WS ]]
            then
                symm_index=$((symm_index+108))
            elif [[ $third == WSe ]]
            then
                symm_index=$((symm_index+144))
            else
                symm_index=$((symm_index+180))
            fi
	    echo "i = "$i" symm_index = "$symm_index
            if (( $i > $symm_index ))
	    then
		echo "already submitted "$M1${X1}2-$M2${X2}2-$M3${X3}2
	    else
		echo "submitting "$M1${X1}2-$M2${X2}2-$M3${X3}2
		submitted=$((submitted+1))
		python snl_prep.py -s $M1${X1}2-$M2${X2}2-$M3${X3}2 -d .
	    fi
            i=$((i+1))
          done
        done
      done
    done
  done
done
echo $submitted
