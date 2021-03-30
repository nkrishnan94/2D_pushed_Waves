

g++ -o outf eff_pop_2d_radial_neutral.cpp -lgsl

for b in 1 2 2.5 3 3.5 4
do
	for u in 75, 125, 150, 175
	do
		for (( i = 0; i < 10 ; i++ )) ### Outer for loop ###
		do
			./outf -B $b -I $u -G $i
		  	cnt=$[cnt+1]
		  	echo $cnt "out of" $tot "\n"

		done
	done
done