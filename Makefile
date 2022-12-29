batchgen:
	python3 batch_gen.py
main: batchgen
	g++ -Wall -pedantic-errors -Wextra -Weffc++ -Wsign-conversion -std=c++17 -O4 -D DEBUG_ASSERTS=0 -D IGNORE_CUDA=1 -D TRACK_TIME=1 -D TRACK_TIME_WALL=1 -g -o main main.cpp -lgmp -lgmpxx
tests: batchgen
	g++ -Wall -pedantic-errors -Wextra -Weffc++ -Wsign-conversion -std=c++17 -D DEBUG_ASSERTS=1 -D IGNORE_CUDA=1 -g -o tests tests.cpp -lgmp -lgmpxx
testbench: batchgen
	# You may want to vary TRACK_TIME
	g++ -Wall -pedantic-errors -Wextra -Weffc++ -Wsign-conversion -std=c++17 -O4 -D DEBUG_ASSERTS=0 -D IGNORE_CUDA=1 -D TRACK_TIME=3 -D TRACK_TIME_WALL=3 -g -o tbm1 testbench_main.cpp
