#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <math.h>
#include <stdbool.h>

#define WIDTH 900
#define HEIGHT 600
#define RAY_NUMBER 400



struct RenderObject 
{
    double x;
    double y;
    int type;
    void (*CollisionDetect)(struct RenderObject*, struct RenderObject*);  // Pointer to collision detection function
    void (*Render)(struct RenderObject*, SDL_Renderer*);
};
struct Circle
{
    struct RenderObject base;
    double radius;
};

//Rays emit at the surface, but from center.
struct LightRay
{
    double x_start,y_start; //position it starts
    double dx, dy; //direction it points. useful for Fill Rays.
};
//Won't implement for now. Would have been nice for color inputs
// void ExtractColors(Uint32 color, Uint8* red, Uint8* green, Uint8* blue, Uint8* alpha){
//     *red = (color >> 24) & 0xFF;  // Extract Red (8 bits)
//     *green = (color >> 16) & 0xFF;  // Extract Green (8 bits)
//     *blue = (color >> 8) & 0xFF;   // Extract Blue (8 bits)
//     *alhpa = color & 0xFF;           // Extract Alpha (8 bits)
// }


void FillCircle(SDL_Renderer* renderer, struct Circle circle){
    double radius_squared = pow(circle.radius,2);
    for (double x = circle.base.x-circle.radius; x<= circle.base.x + circle.radius;x++) {
        for (double y = circle.base.y-circle.radius; y<=circle.base.y +circle.radius;y++){
            double distance_squared = pow(x-circle.base.x,2)+pow(y-circle.base.y,2);
            if (distance_squared  < radius_squared){
                SDL_RenderDrawPoint(renderer, x, y);
            }
        }
    }
}
bool RayCircleIntersection(struct LightRay ray, struct Circle circle, double *t_nearest){
    double dx = ray.dx;
    double dy = ray.dy;
    double x0 = ray.x_start;
    double y0 = ray.y_start;

    double cx = circle.base.x;
    double cy = circle.base.y;
    double r = circle.radius;
    double A = dx * dx + dy * dy;
    double B = 2 * (dx * (x0 - cx) + dy * (y0 - cy));
    double C = (x0 - cx) * (x0 - cx) + (y0 - cy) * (y0 - cy) - r * r;

    double discriminant = B * B - 4 * A * C;

    if (discriminant < 0) {
        return false; // No intersection
    }
    double sqrtD = sqrt(discriminant);
    double t1 = (-B - sqrtD) / (2 * A);
    double t2 = (-B + sqrtD) / (2 * A);

    // Find the nearest valid intersection
    if (t1 >= 0) {
        *t_nearest = t1;
        return true;
    } else if (t2 >= 0) {
        *t_nearest = t2;
        return true;
    }

    return false; // Ray starts inside the circle
}
void FillRays(SDL_Renderer* renderer, struct LightRay rays[RAY_NUMBER], struct Circle occluding_circle){
    for (int i = 0; i< RAY_NUMBER; i++){
        struct LightRay r = rays[i];
        
        double max_expected_distance = 1200; // set to desired circle radius. Set to constant for easy computation
        double x_end = r.x_start + r.dx * max_expected_distance;
        double y_end = r.y_start + r.dy * max_expected_distance;
        double t_nearest;
        if (RayCircleIntersection(r, occluding_circle, &t_nearest)) {
            x_end = r.x_start + r.dx * t_nearest;
            y_end = r.y_start + r.dy * t_nearest;
        }
        SDL_RenderDrawLine(renderer, r.x_start, r.y_start, x_end, y_end);

    }
}

// Stationary Object Collision
// bool CircleCollisionHandle(struct RenderObject* obj1, struct RenderObject* obj2) {
//     struct Circle* circle1 = (struct Circle*) obj1;
//     struct Circle* circle2 = (struct Circle*) obj2;

//     double dx = circle1->base.x - circle2->base.x;
//     double dy = circle1->base.y - circle2->base.y;
    
//     double distance_squared = dx * dx + dy * dy;
//     double radius_sum = circle1->radius + circle2->radius;
//     double radius_sum_squared = radius_sum * radius_sum;

//     // Check if the distance between centers is less than or equal to the sum of radii
//     return distance_squared <= radius_sum_squared;
// }


bool CircleCollisionHandle(struct RenderObject* obj1, struct RenderObject* obj2) {
    struct Circle* circle1 = (struct Circle*)obj1;
    struct Circle* circle2 = (struct Circle*)obj2;

    double dx = circle1->base.x - circle2->base.x;
    double dy = circle1->base.y - circle2->base.y;
    double distance = sqrt(dx * dx + dy * dy);
    double min_distance = circle1->radius + circle2->radius;

    if (distance < min_distance) { // Collision detected
        double overlap = min_distance - distance;  // Amount of overlap

        if (overlap > 0) {
            // Normalize the push direction
            double nx = dx / distance;
            double ny = dy / distance;

            // Push objects apart smoothly by half of the overlap
            circle1->base.x += nx * (overlap / 2.0);
            circle1->base.y += ny * (overlap / 2.0);
            circle2->base.x -= nx * (overlap / 2.0);
            circle2->base.y -= ny * (overlap / 2.0);
        }

        return true;  // Collision happened
    }

    return false;  // No collision
}



bool ObjectCollisionHandle(struct RenderObject* obj1, struct RenderObject* obj2){
    if (obj1->type == 1 && obj2->type ==1) {
        return CircleCollisionHandle(obj1,obj2);
    }
    //add other checks for other objects. 
    return false;
}

void GenerateCircleRays(struct Circle circle, struct LightRay rays[RAY_NUMBER]) {
    // Calculate the angular step size
    double angle_step = 2.0 * M_PI / (double)RAY_NUMBER;
    
    for (int i = 0; i < RAY_NUMBER; i++) {
        // Calculate angle for each ray
        double angle = i * angle_step;

        // Calculate direction vector using cosine and sine
        double dx = cos(angle);
        double dy = sin(angle);

        // Store the ray's start position and direction
        struct LightRay ray = {circle.base.x, circle.base.y, dx, dy};
        rays[i] = ray;

        // For debugging: print out the angle and its direction
        // printf("Ray %d: angle = %f radians (%f degrees), dx = %f, dy = %f\n", i, angle, angle * (180.0 / M_PI), dx, dy);
    }
}


int main() {
    //catch early errors
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        printf("SDL_Init Error: %s\n", SDL_GetError());
        return 1;
    }
    //CREATE WINDOW
    SDL_Window* window = SDL_CreateWindow("Raytracing Demo", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
    if (!window) {
        printf("SDL_CreateWindow Error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    // Create Renderer
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        printf("SDL_CreateRenderer Error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    struct Circle light_circle = {{200, 200, 1},  80}; //Declared circle at (x,y)= (200,200) with r=80

    struct Circle occluding_circle = {{600, 300, 1}, 140};//Declared circle at (x,y)= (500,300) with r=140)
    // animation variable
    struct LightRay rays[RAY_NUMBER];
    // generate initial rays so mouse click doesn't cause them to pop in
    GenerateCircleRays(light_circle,rays);
    int running = 1;
    struct Circle* selected_circle = NULL; // Pointer to track the selected object

    // Event loop to keep the window open. Without this it instantly closes.
    SDL_Event event;
    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT){
                running = 0;
            }
            if (event.type == SDL_MOUSEBUTTONDOWN) {
            double mx = event.button.x;
            double my = event.button.y;

            // Check if the mouse is inside the light circle
            double dx = mx - light_circle.base.x;
            double dy = my - light_circle.base.y;
            if (dx * dx + dy * dy <= light_circle.radius * light_circle.radius) {
                selected_circle = &light_circle;
            }

            // Check if the mouse is inside the occluding circle
            dx = mx - occluding_circle.base.x;
            dy = my - occluding_circle.base.y;
            if (dx * dx + dy * dy <= occluding_circle.radius * occluding_circle.radius) {
                selected_circle = &occluding_circle;
            }
        }

        if (event.type == SDL_MOUSEBUTTONUP) {
            selected_circle = NULL; // Release selection
        }

        if (event.type == SDL_MOUSEMOTION && selected_circle) {
            selected_circle->base.x = event.motion.x;
            selected_circle->base.y = event.motion.y;
            GenerateCircleRays(light_circle, rays); // Recalculate rays for light source
        }
    }

        bool collision = ObjectCollisionHandle((struct RenderObject*)&light_circle, (struct RenderObject*)&occluding_circle);
        if (!collision){
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderClear(renderer);
            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            FillCircle(renderer, light_circle);
            FillCircle(renderer, occluding_circle);
            SDL_SetRenderDrawColor(renderer, 255, 223, 0, 255);//golden yellow
            FillRays(renderer,rays, occluding_circle);
            SDL_RenderPresent(renderer);
        }else{
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);  // Set background to black
            SDL_RenderClear(renderer);  // Clear the screen

            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);  // Set drawing color to white
            FillCircle(renderer, light_circle);  // Draw the main circle
        }
        SDL_Delay(10); // Small delay to reduce CPU usage 100 fps
    }
    //Clear Windows
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
