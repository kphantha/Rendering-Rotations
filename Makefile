all:
	gcc -I src/Include -L src/lib -o rotations rotations.c -lmingw32 -lSDL2main -lSDL2
