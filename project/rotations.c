#include <SDL2/SDL.h>
#include <stdio.h>
#include <math.h>

#define Cw 500
#define Ch 500
#define R1 0.5
#define R2 2
#define ZP 7
#define H 2
#define R3 0.75

void renderFrameTorus(SDL_Renderer *renderer, float A, float B);
void renderFrameArrow(SDL_Renderer *renderer, float A, float B);
void renderFrameCube(SDL_Renderer *renderer, double angle, double q_vec[3]);

int main(int argc, char *argv[])
{
    SDL_Window *window = NULL;
    SDL_Renderer *renderer = NULL;

    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        printf("video initialize error: %s\n", SDL_GetError());
        return 1;
    }
    // set the window and renderer
    if (SDL_CreateWindowAndRenderer(Cw, Ch, 0, &window, &renderer) == -1)
    {
        printf("render initialize error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    double A = 0;
    double B = 0;
    double v[3] = {3, 1, 2};

    for (double frames = 0; frames < 500; frames++)
    {
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        //renderFrameArrow(renderer, A, B);
        renderFrameCube(renderer, A, v);
        SDL_RenderPresent(renderer);
        SDL_Delay(10);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        SDL_PumpEvents();
        A += 0.03;
        // B += 0.01;
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

void renderFrameTorus(SDL_Renderer *renderer, float A, float B)
{
    float distance = Cw * ZP * 3 / (8 * (R1 + R2));

    float cosA = cos(A);
    float sinA = sin(A);
    float cosB = cos(B);
    float sinB = sin(B);

    for (float theta = 0; theta < 2 * M_PI; theta += 0.05)
    {
        float cosTheta = cos(theta);
        float sinTheta = sin(theta);
        // float circlex = R2 + R1 * cosTheta;
        // float circley = R1 * sinTheta;
        for (float phi = 0; phi < 2 * M_PI; phi += 0.05)
        {
            float cosPhi = cos(phi);
            float sinPhi = sin(phi);

            /*Method 1
            double x = circlex * (cosB * cosPhi + sinA * sinB * sinPhi) - circley * cosA * sinB;
            double y = circlex * (sinB * cosPhi - sinA * cosB * sinPhi) + circley * cosA * cosB;
            double z = ZP + cosA * circlex * sinPhi + circley * sinA;
            */

            /* Method 2
            double x = circlex * (cosPhi * cosB - sinB * sinA * sinPhi) + circley * cosA * sinB;
            double y = circlex * (-sinB- * cosPhi - cosB * sinA * sinPhi) + circley * cosA * cosB;
            double z = ZP - circlex * cosA * sinPhi - circley * sinA;
            */

            // Method 3
            double x = R1 * cosTheta * (cosB * cosPhi - sinB * sinA * sinPhi) + R1 * sinB * cosA * sinTheta + R2 * (cosB * cosPhi - sinB * sinA * sinPhi);
            double y = R1 * cosTheta * (-sinB * cosPhi - cosB * sinA * sinPhi) + R1 * cosB * cosA * sinTheta + R2 * (-sinB * cosPhi - cosB * sinA * sinPhi);
            double z = ZP - R2 * cosA * sinPhi - R1 * cosA * sinPhi * cosTheta - R1 * sinA * sinTheta;

            double zInv = 1 / z;
            int xp = (int)(Cw / 2 + distance * zInv * x);
            int yp = (int)(Ch / 2 - distance * zInv * y);

            SDL_RenderDrawPoint(renderer, xp, yp);
        }
    }
}

void renderFrameArrow(SDL_Renderer *renderer, float A, float B)
{
    float distance = Cw * ZP * 2 / (8 * (R1 + R2));

    float cosA = cos(A);
    float sinA = sin(A);
    float cosB = cos(B);
    float sinB = sin(B);

    for (float theta = 0; theta < 2 * M_PI; theta += 0.05)
    {
        float cosTheta = cos(theta);
        float sinTheta = sin(theta);

        for (float height = -H / 2; height < H / 2; height += 0.05)
        {
            double x = R1 * (cosB * cosTheta - sinB * sinA * sinTheta) + height * sinB * cosA;
            double y = R1 * (-sinB * cosTheta - cosB * sinA * sinTheta) + height * cosB * cosA;
            double z = ZP - height * sinA - R1 * cosA * sinTheta;

            double zInv = 1 / z;
            int xp = (int)(Cw / 2 + distance * zInv * x);
            int yp = (int)(Ch / 2 - distance * zInv * y);

            SDL_RenderDrawPoint(renderer, xp, yp);
        }

        float height2 = H / 2;

        for (float radius = R3; radius > 0; radius -= 0.02)
        {
            double x = radius * (cosB * cosTheta - sinB * sinA * sinTheta) + height2 * sinB * cosA;
            double y = radius * (-sinB * cosTheta - cosB * sinA * sinTheta) + height2 * cosB * cosA;
            double z = ZP - height2 * sinA - radius * cosA * sinTheta;

            double zInv = 1 / z;
            int xp = (int)(Cw / 2 + distance * zInv * x);
            int yp = (int)(Ch / 2 - distance * zInv * y);

            SDL_RenderDrawPoint(renderer, xp, yp);

            height2 += 0.03;
        }
    }
}

void renderFrameCube(SDL_Renderer *renderer, double angle, double q_vec[3])
{
    double distance = Cw * ZP * 2 / (8 * H);

    double q_magnitude = sqrt(q_vec[0] * q_vec[0] + q_vec[1] * q_vec[1] + q_vec[2] * q_vec[2]);

    double theta = angle / 2;
    double norm_const = sin(theta) / q_magnitude;

    double a = cos(theta);
    double b = norm_const * q_vec[0];
    double c = norm_const * q_vec[1];
    double d = norm_const * q_vec[2];

    double aa = a * a;
    double bb = b * b;
    double cc = c * c;
    double dd = d * d;
    double ab = 2 * a * b;
    double ac = 2 * a * c;
    double ad = 2 * a * d;
    double bc = 2 * b * c;
    double bd = 2 * b * d;
    double cd = 2 * c * d;

    for (double length = -H / 2; length < H / 2; length += 0.2)
    {
        for (double height = -H / 2; height < H / 2; height += 0.2)
        {
            for (double width = -H / 2; width < H / 2; width += 0.2)
            {
                double x = length * (aa + bb - cc - dd) + height * (bc - ad) + width * (ac + bd);
                double y = length * (ad + bc) + height * (aa - bb + cc - dd) + width * (cd - ab);
                double z = ZP + length * (bd - ac) + height * (ab + cd) + width * (aa - bb - cc + dd);

                double zInv = 1 / z;
                int xp = (int)(Cw / 2 + distance * zInv * x);
                int yp = (int)(Ch / 2 - distance * zInv * y);

                SDL_RenderDrawPoint(renderer, xp, yp);
            }
        }
    }
}