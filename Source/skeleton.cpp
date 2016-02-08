#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* STRUCTS                                                                     */
struct Intersection {
    vec3 position;
    float distance;
    int triangleIndex;
};

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 100;
const int SCREEN_HEIGHT = 100;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;

float focalLength = 50;
vec3 cameraPos(0,0,0);
mat3 R;
float yaw;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle> &triangles, Intersection &closestIntersection);

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    
    LoadTestModel(triangles);

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );

	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

    Uint8* keystate = SDL_GetKeyState(0);
    if(keystate[SDLK_UP]) {
        cameraPos.z += 0.02f;
    } 
    if(keystate[SDLK_DOWN]) {
        cameraPos.z -= 0.02f;
    } 
    if(keystate[SDLK_LEFT]) {
        cameraPos.x += 0.02f;
    } 
    if(keystate[SDLK_RIGHT]) {
        cameraPos.x -= 0.02f;
    } 
}

void Draw()
{

    SDL_FillRect(screen, 0, 0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

    for(int y = 0; y < SCREEN_HEIGHT; y++) {
        for(int x = 0; x < SCREEN_WIDTH; x++) {
            vec3 dir(x-(SCREEN_WIDTH/2), y-(SCREEN_HEIGHT/2),focalLength);
            Intersection closestIntersection;
            if(ClosestIntersection(cameraPos, dir, triangles, closestIntersection)) {
                PutPixelSDL( screen, x, y, triangles[closestIntersection.triangleIndex].color);
            }
            else {
                PutPixelSDL( screen, x, y, vec3(0,0,0));
            }
        }
    }

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle> &triangles, Intersection &closestIntersection) {
    bool foundIntersect = false;

    for(int i=0; i < triangles.size(); i++) {
        Triangle triangle = triangles[i];
        vec3 v0 = triangle.v0;
        vec3 v1 = triangle.v1;
        vec3 v2 = triangle.v2;
        
        vec3 e1 = v1-v0;
        vec3 e2 = v2-v0;
        vec3 b = start-v0;

        mat3 A(-dir,e1,e2);
        vec3 x = glm::inverse(A)*b;

        if(0 < x.y && 0 < x.z && (x.y + x.z) < 1 && 0 <= x.x) {
            float d = glm::distance(x, start);
            if(!foundIntersect || closestIntersection.distance > d) {
                closestIntersection.triangleIndex = i;
                closestIntersection.position = x;
                closestIntersection.distance = d;
                foundIntersect = true;
            }
        }
    }
    return foundIntersect;
}

