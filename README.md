# raytracingC
This application uses SDL2 and C to produce a powerful interactive raytracing demonstration, complete with collision detection, object and wall reflection, distance-based intensity falloff and a toggle for the light rays.

The RayNumber constant can be adjusted to improve performance at the cost of resolved ray number. 

I compiled this project using the following command. Some of the compilation command may depend on your homebrew path variables.

The only required local package installations are SDL2/SDL and the SDL2 TTF module (for text display). You can install them with homebrew using 

brew install sdl2
brew install sdl2_ttf

 gcc raytracing.c -o raytracing \
    -I/opt/homebrew/include \
    -D_REENTRANT \
    -L/opt/homebrew/lib \
    -lSDL2 -lSDL2_ttf

Or single line 

gcc raytracing.c -o raytracing -I/opt/homebrew/include -D_REENTRANT -L/opt/homebrew/lib -lSDL2 -lSDL2_ttf

./raytracing
