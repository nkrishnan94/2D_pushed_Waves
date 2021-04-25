
g++ -o outf eff_pop_2d_radial_neutral.cpp -lgsl

for b in  1.5 2 2.5 3 3.5 4 4.5
do
	for u in 50000 75000 100000 125000 150000 175000 200000 225000
	do
		for (( i = 50; i < 200; i++ )) ### Outer for loop ###
		do
			./outf -B $b -I $u -G $[i*u*2/1000]
		  	cnt=$[cnt+1]
		  	echo $cnt "out of" $tot "\n"

		done
	done

done