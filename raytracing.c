#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <math.h>
#include <stdbool.h>

#define WIDTH 900
#define HEIGHT 600
#define RAY_NUMBER 400
#define OCCLUDINGCIRCLESIZE 80
#define LIGHTCIRCLESIZE 60
#define FALLOFF 2.5
bool light_on = true;


// TODO, add rectangles
struct RenderObject 
{
    double x;
    double y;
    int type;
    bool (*CollisionDetect)(struct RenderObject*, struct RenderObject*);  // Pointer to collision detection function
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
float brightness_buffer[WIDTH][HEIGHT] = {0};  // Zero-initialized

//Won't implement for now. Would have been nice for color inputs
// void ExtractColors(Uint32 color, Uint8* red, Uint8* green, Uint8* blue, Uint8* alpha){
//     *red = (color >> 24) & 0xFF;  // Extract Red (8 bits)
//     *green = (color >> 16) & 0xFF;  // Extract Green (8 bits)
//     *blue = (color >> 8) & 0xFF;   // Extract Blue (8 bits)
//     *alpha = color & 0xFF;           // Extract Alpha (8 bits)
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

void ReflectRayOverCircle(struct LightRay* ray, double hit_x, double hit_y, struct Circle circle) {
    // Compute surface normal at the point of impact
    double nx = hit_x - circle.base.x;
    double ny = hit_y - circle.base.y;

    // Normalize the normal
    double length = sqrt(nx * nx + ny * ny);
    if (length == 0) return; // Avoid divide by zero
    nx /= length;
    ny /= length;

    // Compute dot product of ray direction and normal
    double dot = ray->dx * nx + ray->dy * ny;

    // Reflect direction over normal
    ray->dx = ray->dx - 2 * dot * nx;
    ray->dy = ray->dy - 2 * dot * ny;

    // Normalize the new direction
    length = sqrt(ray->dx * ray->dx + ray->dy * ray->dy);
    ray->dx /= length;
    ray->dy /= length;
}

void FillBouncingRays(SDL_Renderer* renderer, struct LightRay rays[RAY_NUMBER], struct Circle occluder, struct Circle light_source, double falloff_factor) {

    const double max_distance = 1200.0;
    const double intensity_threshold = 0.01;
    const int max_bounces = 10;

    for (int i = 0; i < RAY_NUMBER; i++) {
        double x = rays[i].x_start;
        double y = rays[i].y_start;
        double dx = rays[i].dx;
        double dy = rays[i].dy;

        double traveled = 0.0;
        int bounces = 0;

        while (bounces < max_bounces) {
            // Step 1: Determine distance to each wall
            double tx = dx > 0 ? (WIDTH - x) / dx : (dx < 0 ? (0 - x) / dx : 1e9);
            double ty = dy > 0 ? (HEIGHT - y) / dy : (dy < 0 ? (0 - y) / dy : 1e9);
            double t_wall = fmin(tx, ty);
            // Step 2: Determine distance to occluders
            double t_circle = 1e9;
            double t_nearest = 1e9;
            bool hit_light = false;

            struct LightRay ray = {x, y, dx, dy};

            // Check occluding object
            if (RayCircleIntersection(ray, occluder, &t_nearest) && t_nearest > 0) {
                t_circle = t_nearest;
            }

            // Check light source if not the first segment
            if (bounces > 0) {
                double t_light;
                if (RayCircleIntersection(ray, light_source, &t_light) && t_light > 0) {
                    if (t_light < t_circle) {
                        t_circle = t_light;
                        hit_light = true;  // Mark that we hit the light source
                    }
                }
            }





            // Step 3: Choose minimum valid t
            double t_min = fmin(t_wall, t_circle);

            // Step 4: Clip to remaining max distance
            double segment_length = t_min;
            if (traveled + segment_length > max_distance) {
                segment_length = max_distance - traveled;
            }

            double x2 = x + dx * segment_length;
            double y2 = y + dy * segment_length;

            //Makes a weird light source shadow effect
            // // If the ray hits the light source on a later bounce, terminate
            // if (bounces > 0 && hit_light) {
            //     break;
            // }

            // Step 5: Compute falloff at segment start and end
            double t0 = traveled / max_distance;
            double t1 = (traveled + segment_length) / max_distance;
            double intensity0 = exp(-falloff_factor * t0);
            double intensity1 = exp(-falloff_factor * t1);

            if (intensity0 < intensity_threshold) break;

            Uint8 r0 = (Uint8)(255 * intensity0);
            Uint8 g0 = (Uint8)(223 * intensity0);
            Uint8 r1 = (Uint8)(255 * intensity1);
            Uint8 g1 = (Uint8)(223 * intensity1);

            int steps = (int)(segment_length / 1.5);
            if (steps < 1) steps = 1;

            for (int s = 0; s < steps; s++) {
                double s0 = (double)s / steps;
                double s1 = (double)(s + 1) / steps;

                double xs = x + dx * segment_length * s0;
                double ys = y + dy * segment_length * s0;
                double xe = x + dx * segment_length * s1;
                double ye = y + dy * segment_length * s1;

                double dist_start = traveled + segment_length * s0;
                double dist_end   = traveled + segment_length * s1;

                double intensity_start = exp(-falloff_factor * (dist_start / max_distance));
                double intensity_end   = exp(-falloff_factor * (dist_end   / max_distance));
                if (intensity_start < intensity_threshold) break;

                Uint8 r = (Uint8)(255 * (intensity_start + intensity_end) / 2);
                Uint8 g = (Uint8)(223 * (intensity_start + intensity_end) / 2);

                int ix0 = (int)xs;
                int iy0 = (int)ys;
                int ix1 = (int)xe;
                int iy1 = (int)ye;

                if (ix0 >= 0 && ix0 < WIDTH && iy0 >= 0 && iy0 < HEIGHT &&
                    ix1 >= 0 && ix1 < WIDTH && iy1 >= 0 && iy1 < HEIGHT) {
                    float avg_intensity = (intensity_start + intensity_end) / 2.0f;


                    bool can_draw = true;

                    // Check all pixels on the line before drawing. a little slow but massively improved visual fidelity.
                    int dx = abs(ix1 - ix0);
                    int dy = abs(iy1 - iy0);
                    int sx = ix0 < ix1 ? 1 : -1;
                    int sy = iy0 < iy1 ? 1 : -1;
                    int err = dx - dy;
                    int x = ix0, y = iy0;

                    while (true) {
                        if (x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT) {
                            if (avg_intensity <= brightness_buffer[x][y]) {
                                can_draw = false;
                                break;
                            }
                        }
                        if (x == ix1 && y == iy1) break;
                        int e2 = 2 * err;
                        if (e2 > -dy) { err -= dy; x += sx; }
                        if (e2 < dx)  { err += dx; y += sy; }
                    }

                    if (can_draw) {
                        // Commit pixel updates and draw line
                        x = ix0; y = iy0; err = dx - dy;
                        while (true) {
                            if (x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT)
                                brightness_buffer[x][y] = avg_intensity;
                            if (x == ix1 && y == iy1) break;
                            int e2 = 2 * err;
                            if (e2 > -dy) { err -= dy; x += sx; }
                            if (e2 < dx)  { err += dx; y += sy; }
                        }

                        SDL_SetRenderDrawColor(renderer, r, g, 0, 255);
                        SDL_RenderDrawLine(renderer, ix0, iy0, ix1, iy1);
                    }

                }


            }


            // Step 6: Update state for next segment or exit
            traveled += segment_length;
            x = x2;
            y = y2;

            // if (t_circle <= t_wall) break; // HIT THE OCCLUDING CIRCLE — STOP
            if (t_circle <= t_wall){
                ReflectRayOverCircle(&ray, x2, y2, occluder);
                dx = ray.dx;
                dy = ray.dy;
                // Stop if ray hit the occluder and we're not on the first bounce
                if (bounces > 0 && t_circle <= t_wall) {
                    break;
                }
            } else {
                // reflect off walls
                if (t_wall == tx) dx = -dx;
                if (t_wall == ty) dy = -dy;
            }
            



            bounces++;
        }
    }
}



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

        // printf("Ray %d: angle = %f radians (%f degrees), dx = %f, dy = %f\n", i, angle, angle * (180.0 / M_PI), dx, dy);
    }
}


int main() {
    //catch early errors
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        printf("SDL_Init Error: %s\n", SDL_GetError());
        return 1;
    }
    //initialize ttf 
    if (TTF_Init() != 0) {
    printf("TTF_Init Error: %s\n", TTF_GetError());
    SDL_Quit();
    return 1;
    }

    //CREATE WINDOW
    SDL_Window* window = SDL_CreateWindow("Raytracing Demo", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
    if (!window) {
        printf("SDL_CreateWindow Error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    //Light Switch Button
    SDL_Rect button_rect = {10, 10, 120, 40};

    // Create Renderer
    // SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED );
    if (!renderer) {
        printf("SDL_CreateRenderer Error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    TTF_Font* font = TTF_OpenFont("/Library/Fonts/Arial Unicode.ttf", 20); // or your system's path to a .ttf
    if (!font) {
        printf("TTF_OpenFont Error: %s\n", TTF_GetError());
        SDL_Quit();
        return 1;
    }

    struct Circle light_circle = {{200, 200, 1},  LIGHTCIRCLESIZE}; //Declared circle at (x,y)= (200,200) with r=80
    light_circle.base.CollisionDetect = ObjectCollisionHandle;
    
    
    struct Circle occluding_circle = {{600, 300, 1}, OCCLUDINGCIRCLESIZE};
    occluding_circle.base.CollisionDetect = ObjectCollisionHandle;
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
                if (event.button.button == SDL_BUTTON_LEFT) {
                    int mx = event.button.x;
                    int my = event.button.y;

                    if (mx >= button_rect.x && mx <= button_rect.x + button_rect.w &&
                        my >= button_rect.y && my <= button_rect.y + button_rect.h) {
                        light_on = !light_on;
                        if (light_on) {
                            GenerateCircleRays(light_circle, rays);
                        } else {
                            // Set all rays to zero vectors
                            for (int i = 0; i < RAY_NUMBER; i++) {
                                rays[i] = (struct LightRay){0, 0, 0, 0};
                            }
                        }
                    }
                }
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
            double new_x = event.motion.x;
            double new_y = event.motion.y;

            // Clamp to window bounds
            double r = selected_circle->radius;
            if (new_x < r) new_x = r;
            if (new_x > WIDTH - r) new_x = WIDTH - r;
            if (new_y < r) new_y = r;
            if (new_y > HEIGHT - r) new_y = HEIGHT - r;


            selected_circle->base.x = new_x;
            selected_circle->base.y = new_y;

            if (selected_circle == &light_circle)
                GenerateCircleRays(light_circle, rays); // Recalculate rays only if light moved
        }

    }

        // bool collision = ObjectCollisionHandle((struct RenderObject*)&light_circle, (struct RenderObject*)&occluding_circle);
        // bool collision = false;
        bool collision = light_circle.base.CollisionDetect((struct RenderObject*)&light_circle, (struct RenderObject*)&occluding_circle);
        // clear and draw scene
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        for (int x = 0; x < WIDTH; x++) {
            for (int y = 0; y < HEIGHT; y++) {
                brightness_buffer[x][y] = 0.0f;
            }
        }
        SDL_RenderClear(renderer);

        // draw objects
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        FillCircle(renderer, light_circle);
        FillCircle(renderer, occluding_circle);

        // draw rays only if light is on (don't care about collision status)
        if (light_on) {
            SDL_SetRenderDrawColor(renderer, 255, 223, 0, 255);
            FillBouncingRays(renderer, rays, occluding_circle, light_circle, FALLOFF);
        }

        // Always draw the toggle button
        SDL_SetRenderDrawColor(renderer, 100, 100, 100, 255);
        SDL_RenderFillRect(renderer, &button_rect);

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderDrawRect(renderer, &button_rect);
        // Create surface and texture for label
        SDL_Color textColor = {255, 255, 255};
        const char* label_text = light_on ? "LIGHT ON" : "LIGHT OFF";

        SDL_Surface* textSurface = TTF_RenderText_Blended(font, label_text, textColor);
        if (textSurface) {
            SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
            if (textTexture) {
                int textW, textH;
                SDL_QueryTexture(textTexture, NULL, NULL, &textW, &textH);
                SDL_Rect textRect = {
                    button_rect.x + (button_rect.w - textW) / 2,
                    button_rect.y + (button_rect.h - textH) / 2,
                    textW,
                    textH
                };
                SDL_RenderCopy(renderer, textTexture, NULL, &textRect);
                SDL_DestroyTexture(textTexture);
            }
            SDL_FreeSurface(textSurface);
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(10); // Small delay to reduce CPU usage 100 fps
    }
    //Clear Windows
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
