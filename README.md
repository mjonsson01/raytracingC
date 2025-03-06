# raytracingC
This application uses SDL2 and C to produce a rudimentary interactive raytracing demonstration, complete with collision detection.

The RayNumber constant can be adjusted to improve performance at the cost of resolved ray number. The collision detection function CircleCollisionHandle can be replaced with the commmented-out function to create a stationary contact object.

I compiled this project using. Some of the compilation command will depend on your homebrew path variables.

The only required local package installations are SDL2/SDL. You can install them with homebrew using 

brew install sdl2

gcc raytracing.c -o raytracing -I/opt/homebrew/include -D_REENTRANT -L/opt/homebrew/lib -lSDL2

./raytracing
