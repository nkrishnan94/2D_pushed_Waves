
g++ -o outf eff_pop_2d_radial_neutral.cpp -lgsl

for b in  0 2.5 3 3.5 5  
do
	for u in 30 40 50 60 70 80 90
	do
		for (( i = 100; i < 200; i++ )) ### Outer for loop ###
		do
			./outf -B $b -U $u -G $i
		  	cnt=$[cnt+1]
		  	echo $cnt "out of" $tot "\n"

		done
	done

done