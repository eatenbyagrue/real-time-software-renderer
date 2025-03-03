#include <SDL.h>
#include <stdio.h>

void loop() {
    SDL_Event e;
    while (SDL_PollEvent(&e) != 0) {
        if (e.type == SDL_KEYDOWN) {
            if (e.key.keysym.sym == SDLK_a) {
                printf("When I press the A key on my keyboard, this branch does NOT execute\n");
            } else if (e.key.keysym.sym == SDL_SCANCODE_A) {
                printf("This branch DOES execute when I press the A key on my keyboard.\n");
            }

            printf("sym %i scancode %i\n", e.key.keysym.sym, e.key.keysym.scancode);
        }
    }
}

int main(int argc, char **argv) {
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *window =
        SDL_CreateWindow("sdl2_key", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 640, 480, 0);

    while (true) {
        loop();
    }

    return 0;
}
