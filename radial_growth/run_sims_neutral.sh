
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=125




for b in 1 2 3 5 7
do
      for (( i = 0; i < 25 ; i++ )) ### Outer for loop ###
      do
      	g++ -o outf eff_pop_2d_radial_neutral.cpp -lgsl
      	./outf -B $b -U 100
      	cnt=$[cnt+1]
      	echo $cnt "out of" $tot "\n"
      done 	
   done
done



