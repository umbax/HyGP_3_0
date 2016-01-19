# name of the executable
OUT = z

${OUT} : TINL2_mod.o MINL2.o
	 g++ TINL2_mod.o MINL2.o -L./ -lf2c -lm -o ${OUT} #linking

TINL2.o : TINL2_mod.cpp 
	g++ -c TINL2_mod.cpp -o TINL2_mod.o  # -c compile without linking

MINL2.o : MINL2.cpp
	g++ -c MINL2.cpp -o MINL2.o # -c compile without linking

clean:	#this line should contain dependencies. As there are not, leave it blank!
				echo "Cleaning up files"
				rm TINL2_mod.o
				rm MINL2.o
				rm ${OUT}
				echo "Cleaning done"
